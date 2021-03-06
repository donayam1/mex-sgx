/*
 * Copyright (C) 2011-2021 Intel Corporation. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of Intel Corporation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

# include <unistd.h>
# include <pwd.h>
# define MAX_PATH FILENAME_MAX

#include <sgx_urts.h>
#include <sgx_uswitchless.h>
#include "App.h"
#include "Enclave_u.h"

/////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <limits>
#include <unistd.h>
#include <omp.h>
#include <fstream>
#include "utils.h"
#include "../Aes/aes.hpp"

// #define VTUNE_ANALYSIS 1

#define CLMUL 8
#ifdef VTUNE_ANALYSIS
    #include <ittnotify.h>
#endif

#define DEFAULT_MATCH 1
#define DEFAULT_MISMATCH 4
#define DEFAULT_OPEN 6
#define DEFAULT_EXTEND 1
#define DEFAULT_AMBIG -1

#undef MAX_SEQ_LEN_REF
#define MAX_SEQ_LEN_REF 2048
#undef MAX_SEQ_LEN_QER
#define MAX_SEQ_LEN_QER 256

///////////////////////////////////////////////////////////////
#if defined(__clang__) || defined(__GNUC__)
#define __mmask8 uint8_t
#define __mmask16 uint16_t
#define __mmask32 uint32_t
#endif

// ------------------------------------------------------------------------------------
// MACROs for vector code
extern uint64_t prof[10][112];
#define AMBIG 4
#define DUMMY1 99
#define DUMMY2 100

///////////////////////////////////////////////////////////////////////////////////////////


/* Global EID shared by multiple threads */
sgx_enclave_id_t global_eid = 0;

typedef struct _sgx_errlist_t {
    sgx_status_t err;
    const char *msg;
    const char *sug; /* Suggestion */
} sgx_errlist_t;

#define REPEATS 500000

/* Error code returned by sgx_create_enclave */
static sgx_errlist_t sgx_errlist[] = {
    {
        SGX_ERROR_UNEXPECTED,
        "Unexpected error occurred.",
        NULL
    },
    {
        SGX_ERROR_INVALID_PARAMETER,
        "Invalid parameter.",
        NULL
    },
    {
        SGX_ERROR_OUT_OF_MEMORY,
        "Out of memory.",
        NULL
    },
    {
        SGX_ERROR_ENCLAVE_LOST,
        "Power transition occurred.",
        "Please refer to the sample \"PowerTransition\" for details."
    },
    {
        SGX_ERROR_INVALID_ENCLAVE,
        "Invalid enclave image.",
        NULL
    },
    {
        SGX_ERROR_INVALID_ENCLAVE_ID,
        "Invalid enclave identification.",
        NULL
    },
    {
        SGX_ERROR_INVALID_SIGNATURE,
        "Invalid enclave signature.",
        NULL
    },
    {
        SGX_ERROR_OUT_OF_EPC,
        "Out of EPC memory.",
        NULL
    },
    {
        SGX_ERROR_NO_DEVICE,
        "Invalid SGX device.",
        "Please make sure SGX module is enabled in the BIOS, and install SGX driver afterwards."
    },
    {
        SGX_ERROR_MEMORY_MAP_CONFLICT,
        "Memory map conflicted.",
        NULL
    },
    {
        SGX_ERROR_INVALID_METADATA,
        "Invalid enclave metadata.",
        NULL
    },
    {
        SGX_ERROR_DEVICE_BUSY,
        "SGX device was busy.",
        NULL
    },
    {
        SGX_ERROR_INVALID_VERSION,
        "Enclave version was invalid.",
        NULL
    },
    {
        SGX_ERROR_INVALID_ATTRIBUTE,
        "Enclave was not authorized.",
        NULL
    },
    {
        SGX_ERROR_ENCLAVE_FILE_ACCESS,
        "Can't open enclave file.",
        NULL
    },
};

///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////
/* Check error conditions for loading enclave */
void print_error_message(sgx_status_t ret)
{
    size_t idx = 0;
    size_t ttl = sizeof sgx_errlist/sizeof sgx_errlist[0];

    for (idx = 0; idx < ttl; idx++) {
        if(ret == sgx_errlist[idx].err) {
            if(NULL != sgx_errlist[idx].sug)
                printf("Info: %s\n", sgx_errlist[idx].sug);
            printf("Error: %s\n", sgx_errlist[idx].msg);
            break;
        }
    }

    if (idx == ttl)
        printf("Error: Unexpected error occurred.\n");
}

/* Initialize the enclave:
 *   Call sgx_create_enclave to initialize an enclave instance
 */
int initialize_enclave(const sgx_uswitchless_config_t* us_config)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;

    /* Call sgx_create_enclave to initialize an enclave instance */
    /* Debug Support: set 2nd parameter to 1 */

    const void* enclave_ex_p[32] = { 0 };

    enclave_ex_p[SGX_CREATE_ENCLAVE_EX_SWITCHLESS_BIT_IDX] = (const void*)us_config;

    ret = sgx_create_enclave_ex(ENCLAVE_FILENAME, SGX_DEBUG_FLAG, NULL, NULL, &global_eid, NULL, SGX_CREATE_ENCLAVE_EX_SWITCHLESS, enclave_ex_p);
    if (ret != SGX_SUCCESS) {
        print_error_message(ret);
        return -1;
    }

    return 0;
}

/* OCall functions */
void ocall_print_string(const char *str,size_t length)
{
    /* Proxy/Bridge will check the length and null-terminate
     * the input string to prevent buffer overflow.
     */
    printf("%s\n", str);
}

// void ocall_return(int maa) {
//     printf("\n\n The retuned value is %d \n\n\n", maa);
// }

//void ocall_empty(void) {}
//void ocall_empty_switchless(void) {}

/* Application entry */


/*
uint8_t data[]={0x0,0x46,0xb5,0x28,0x58,0x3c,0x91,0x4b,0xd2,0x5c,0x25,0x80,0x30,0xce,0xd3,0x8f,'\0'};
ATTCGATTCCAGCGTA

0xf6,0x49,0xf6,0x23,0xe8,0x14,0xba,0x55,0x8e,0x1f,0x58,0x75,0xbb,0x9f,0x73,0x90
GCTCGATCCCAGCGTA

0x9a,0xc1,0xd5,0xb8,0x39,0x93,0x4c,0xf0,0xd9,0x25,0x45,0xe4,0xa,0xeb,0x32,0x1c
GCATGATCCCATGGTC
*/

int SGX_CDECL main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;

// Configuration for Switchless SGX 
    sgx_uswitchless_config_t us_config = SGX_USWITCHLESS_CONFIG_INITIALIZER;
    us_config.num_uworkers = 2;
    us_config.num_tworkers = 2;

 // Initialize the enclave 
    if(initialize_enclave(&us_config) < 0)
    {
        printf("Error: enclave initialization failed\n");
        return -1;
    }

    // char data[32] = "ACTGACTGACTGACT";
    // int length = strlen(data);
    // printf("Length = %d \n",length);
    // while(length <=30){
    //     strncat(data,"0",1);
    //     length = strlen(data);
    //     printf("Length %d= %s \n",length,data);
       
    // }

    // struct AES_ctx ctx;
    // AES_init_ctx(&ctx, key);

    uint8_t data[]={0xa0,0x26,0xc0,0x5f,0x26,0x84,0x5c,0x56,0xce,0x20,0x2,0xe8,0x36,0xe8,0x40,0x1f,'\0'};
    
    // AES_ECB_decrypt(&ctx, data);
    // char *datac = (char*)data;
    // printf ("App decrypted = %s",datac);

    const char *input = (char*)data;//argv[1]; 
    size_t size = 17;//strlen(input);//argv[1]
    printf("size=%d\n",size);
    int matchCount;
    sgx_status_t status1 = ecall_matchOligo(global_eid,&matchCount,data,size);
    if (status1 != SGX_SUCCESS) {
        printf("ERROR: ECall failed\n");
        print_error_message(status1);
        exit(-1);
    }

    printf("Done. with match count %d\n",matchCount);



    sgx_destroy_enclave(global_eid);

//////////////////////////////////////////////////////////////////////////

#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
	
    return 0;
}


