
#include "mex.h"
#include "../Enclave_t.h"
#include "ltable.h"
#include <string.h>
// #include "input.h"

// #include "server.h"

// LookupTable::LookupTable()
// {
// }
// LookupTable::LookupTable() //SequenceSet *g
// {
//   word_length = WORD_LENGTH;
//   // genome = g;
//   // Default command-line parameter values
//   mismatches = 1;
//   mismatches2 = 2;
//   indels = 1;
//   server_port = PORT;
//   genome_size = MAX_SEQUENCE_LENGTH;
//   // strcpy(xfilein, "");
//   // strcpy(xfileout, "");
//   ltable_step = 1;
//   max_hits = 100;
// }
  // char oligo[MAX_CHAR];
  int word_length = WORD_LENGTH,ltable_step=1;
  
  // int ltable_count[MAX_LTABLE];
  // int *ltable[MAX_LTABLE];

  int global_pos;
  int c;
  char match_sense[MAX_MATCHES];
  int match_shift[MAX_MATCHES], matches[MAX_MATCHES];
  int query_len, value, pos;
  int query_full_len;
  int sense_change,match_count; //j
  // Parameters set from command lines
  int mismatches=1, mismatches2=2, indels=1, max_hits=100;//,  genome_size=MAX_SEQUENCE_LENGTH;
  char primer_label[128];//,xfilein[128],xfileout[128];




  char visited[MAX_LTABLE], word[MAX_CHAR], new_word[MAX_CHAR];
  int Power_of_4[12] = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576};
  //  int mismatch_depth[3] = {1,34,530};
  int mismatch_depth[3] = {1, 34, 528};
  int mismatch_base[532] = {
      0, 1, 2, 3, 4, 8, 12, 16, 32, 48, 64, 128, 192, 256, 512, 768, 1024,
      2048, 3072, 4096, 8192, 12288, 16384, 32768, 49152, 65536,
      131072, 196608, 262144, 524288, 786432, 1048576, 2097152, 3145728,
      5, 9, 13, 17, 33, 49, 65, 129, 193, 257, 513, 769, 1025,
      2049, 3073, 4097, 8193, 12289, 16385, 32769, 49153, 65537,
      131073, 196609, 262145, 524289, 786433, 1048577, 2097153, 3145729,
      6, 10, 14, 18, 34, 50, 66, 130, 194, 258, 514, 770, 1026,
      2050, 3074, 4098, 8194, 12290, 16386, 32770, 49154, 65538,
      131074, 196610, 262146, 524290, 786434, 1048578, 2097154, 3145730,
      7, 11, 15, 19, 35, 51, 67, 131, 195, 259, 515, 771, 1027,
      2051, 3075, 4099, 8195, 12291, 16387, 32771, 49155, 65539,
      131075, 196611, 262147, 524291, 786435, 1048579, 2097155, 3145731,
      20, 36, 52, 68, 132, 196, 260, 516, 772, 1028,
      2052, 3076, 4100, 8196, 12292, 16388, 32772, 49156, 65540,
      131076, 196612, 262148, 524292, 786436, 1048580, 2097156, 3145732,
      24, 40, 56, 72, 136, 200, 264, 520, 776, 1032,
      2056, 3080, 4104, 8200, 12296, 16392, 32776, 49160, 65544,
      131080, 196616, 262152, 524296, 786440, 1048584, 2097160, 3145736,
      28, 44, 60, 76, 140, 204, 268, 524, 780, 1036,
      2060, 3084, 4108, 8204, 12300, 16396, 32780, 49164, 65548,
      131084, 196620, 262156, 524300, 786444, 1048588, 2097164, 3145740,
      80, 144, 208, 272, 528, 784, 1040,
      2064, 3088, 4112, 8208, 12304, 16400, 32784, 49168, 65552,
      131088, 196624, 262160, 524304, 786448, 1048592, 2097168, 3145744,
      96, 160, 224, 288, 544, 800, 1056,
      2080, 3104, 4128, 8224, 12320, 16416, 32800, 49184, 65568,
      131104, 196640, 262176, 524320, 786464, 1048608, 2097184, 3145760,
      112, 176, 240, 304, 560, 816, 1072,
      2096, 3120, 4144, 8240, 12336, 16432, 32816, 49200, 65584,
      131120, 196656, 262192, 524336, 786480, 1048624, 2097200, 3145776,
      320, 576, 832, 1088,
      2112, 3136, 4160, 8256, 12352, 16448, 32832, 49216, 65600,
      131136, 196672, 262208, 524352, 786496, 1048640, 2097216, 3145792,
      384, 640, 896, 1152,
      2176, 3200, 4224, 8320, 12416, 16512, 32896, 49280, 65664,
      131200, 196736, 262272, 524416, 786560, 1048704, 2097280, 3145856,
      448, 704, 960, 1216,
      2240, 3264, 4288, 8384, 12480, 16576, 32960, 49344, 65728,
      131264, 196800, 262336, 524480, 786624, 1048768, 2097344, 3145920,
      1280, 2304, 3328, 4352, 8448, 12544, 16640, 33024, 49408, 65792,
      131328, 196864, 262400, 524544, 786688, 1048832, 2097408, 3145984,
      1536, 2560, 3584, 4608, 8704, 12800, 16896, 33280, 49664, 66048,
      131584, 197120, 262656, 524796, 786944, 1049088, 2097664, 3146240,
      1792, 2816, 3840, 4864, 8960, 13056, 17152, 33536, 49920, 66304,
      131840, 197376, 262912, 525056, 787200, 1049344, 2097920, 3146496,
      5120, 9216, 13312, 17408, 33792, 50176, 66560,
      132096, 197632, 263168, 525312, 787456, 1049600, 2098176, 3146752,
      6144, 10240, 14336, 18432, 34816, 51200, 67584,
      133120, 198656, 264192, 526336, 788480, 1050624, 2099200, 3147776,
      7168, 11264, 15360, 19456, 35840, 52224, 68608,
      134144, 199680, 265216, 527360, 789504, 1051648, 2100224, 3148800,
      20480, 36864, 53248, 69632,
      135168, 200704, 266240, 528384, 790528, 1052672, 2101248, 3149824,
      24576, 40960, 57344, 73728,
      139264, 204800, 270336, 532480, 794624, 1056768, 2105344, 3153920,
      28672, 45056, 61440, 77824,
      143360, 208896, 274432, 536576, 798720, 1060864, 2109440, 3158016,
      81920, 147456, 212992, 278528, 540672, 802816, 1064960, 2113536, 3162112,
      98304, 163840, 229376, 294912, 557056, 819200, 1081344, 2129920, 3178496,
      114688, 180224, 245760, 311296, 573440, 835584, 1097728, 2146304, 3194880,
      327680, 589824, 851968, 1114112, 2162688, 3211264,
      393216, 655360, 917504, 1179648, 2228224, 3276800,
      458752, 720896, 983040, 1245184, 2293760, 3342336,
      1310720, 2359296, 3407872,
      1572864, 2621440, 3670016,
      1835008, 2883584, 3932160};
  int mask1[10] = {
      1048575, 1048572, 1048560, 1048512, 1048320, 1047552, 1044480, 1032192,
      983040, 786432};
  int mask2[10] = {
      3, 15, 63, 255, 1023, 4095, 16383, 65535, 262143, 1048575};






int ecall_matchOligo(char *oligo,size_t size)
{
  
  query_full_len = strlen(oligo);
  ocall_print_string(oligo,query_full_len);
  // char *nl="\n";
  // ocall_print_string(nl,strlen(nl));
  // justCall(oligo);

  find_seed_matches(oligo);
  sense_change = match_count;
  //cout << "Executing find_matches reverse." << endl;
  char reversed_string[MAX_CHAR];
  reverse(oligo,reversed_string,MAX_CHAR);
  query_full_len = strlen(reversed_string);
  find_seed_matches(reversed_string);
  int j;
  for (j = 0; j < sense_change; j++)
  {
    match_sense[j] = '+';
  }
  for (j = sense_change; j < match_count; j++)
  {
    match_sense[j] = '-';
  }
  return match_count;
}
// void justCall(char *string){
//   char *st = "another called\n";
//   ocall_print_string(st,strlen(st));
// }
void find_seed_matches(char *string)
// void find_seed_matches(char *string, char *result)
{
  //cout << "Entering find_matches." << endl;
  // char *st = "another called\n";
  // ocall_print_string(st,strlen(st));

  int bin, bin1, bin2;
  int q, base, segment_step, shift, ltable_shift, lts;
  // char visited[MAX_LTABLE], word[MAX_CHAR], new_word[MAX_CHAR];
  // int Power_of_4[12] = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576};
  // //  int mismatch_depth[3] = {1,34,530};
  // int mismatch_depth[3] = {1, 34, 528};
  // int mismatch_base[532] = {
  //     0, 1, 2, 3, 4, 8, 12, 16, 32, 48, 64, 128, 192, 256, 512, 768, 1024,
  //     2048, 3072, 4096, 8192, 12288, 16384, 32768, 49152, 65536,
  //     131072, 196608, 262144, 524288, 786432, 1048576, 2097152, 3145728,
  //     5, 9, 13, 17, 33, 49, 65, 129, 193, 257, 513, 769, 1025,
  //     2049, 3073, 4097, 8193, 12289, 16385, 32769, 49153, 65537,
  //     131073, 196609, 262145, 524289, 786433, 1048577, 2097153, 3145729,
  //     6, 10, 14, 18, 34, 50, 66, 130, 194, 258, 514, 770, 1026,
  //     2050, 3074, 4098, 8194, 12290, 16386, 32770, 49154, 65538,
  //     131074, 196610, 262146, 524290, 786434, 1048578, 2097154, 3145730,
  //     7, 11, 15, 19, 35, 51, 67, 131, 195, 259, 515, 771, 1027,
  //     2051, 3075, 4099, 8195, 12291, 16387, 32771, 49155, 65539,
  //     131075, 196611, 262147, 524291, 786435, 1048579, 2097155, 3145731,
  //     20, 36, 52, 68, 132, 196, 260, 516, 772, 1028,
  //     2052, 3076, 4100, 8196, 12292, 16388, 32772, 49156, 65540,
  //     131076, 196612, 262148, 524292, 786436, 1048580, 2097156, 3145732,
  //     24, 40, 56, 72, 136, 200, 264, 520, 776, 1032,
  //     2056, 3080, 4104, 8200, 12296, 16392, 32776, 49160, 65544,
  //     131080, 196616, 262152, 524296, 786440, 1048584, 2097160, 3145736,
  //     28, 44, 60, 76, 140, 204, 268, 524, 780, 1036,
  //     2060, 3084, 4108, 8204, 12300, 16396, 32780, 49164, 65548,
  //     131084, 196620, 262156, 524300, 786444, 1048588, 2097164, 3145740,
  //     80, 144, 208, 272, 528, 784, 1040,
  //     2064, 3088, 4112, 8208, 12304, 16400, 32784, 49168, 65552,
  //     131088, 196624, 262160, 524304, 786448, 1048592, 2097168, 3145744,
  //     96, 160, 224, 288, 544, 800, 1056,
  //     2080, 3104, 4128, 8224, 12320, 16416, 32800, 49184, 65568,
  //     131104, 196640, 262176, 524320, 786464, 1048608, 2097184, 3145760,
  //     112, 176, 240, 304, 560, 816, 1072,
  //     2096, 3120, 4144, 8240, 12336, 16432, 32816, 49200, 65584,
  //     131120, 196656, 262192, 524336, 786480, 1048624, 2097200, 3145776,
  //     320, 576, 832, 1088,
  //     2112, 3136, 4160, 8256, 12352, 16448, 32832, 49216, 65600,
  //     131136, 196672, 262208, 524352, 786496, 1048640, 2097216, 3145792,
  //     384, 640, 896, 1152,
  //     2176, 3200, 4224, 8320, 12416, 16512, 32896, 49280, 65664,
  //     131200, 196736, 262272, 524416, 786560, 1048704, 2097280, 3145856,
  //     448, 704, 960, 1216,
  //     2240, 3264, 4288, 8384, 12480, 16576, 32960, 49344, 65728,
  //     131264, 196800, 262336, 524480, 786624, 1048768, 2097344, 3145920,
  //     1280, 2304, 3328, 4352, 8448, 12544, 16640, 33024, 49408, 65792,
  //     131328, 196864, 262400, 524544, 786688, 1048832, 2097408, 3145984,
  //     1536, 2560, 3584, 4608, 8704, 12800, 16896, 33280, 49664, 66048,
  //     131584, 197120, 262656, 524796, 786944, 1049088, 2097664, 3146240,
  //     1792, 2816, 3840, 4864, 8960, 13056, 17152, 33536, 49920, 66304,
  //     131840, 197376, 262912, 525056, 787200, 1049344, 2097920, 3146496,
  //     5120, 9216, 13312, 17408, 33792, 50176, 66560,
  //     132096, 197632, 263168, 525312, 787456, 1049600, 2098176, 3146752,
  //     6144, 10240, 14336, 18432, 34816, 51200, 67584,
  //     133120, 198656, 264192, 526336, 788480, 1050624, 2099200, 3147776,
  //     7168, 11264, 15360, 19456, 35840, 52224, 68608,
  //     134144, 199680, 265216, 527360, 789504, 1051648, 2100224, 3148800,
  //     20480, 36864, 53248, 69632,
  //     135168, 200704, 266240, 528384, 790528, 1052672, 2101248, 3149824,
  //     24576, 40960, 57344, 73728,
  //     139264, 204800, 270336, 532480, 794624, 1056768, 2105344, 3153920,
  //     28672, 45056, 61440, 77824,
  //     143360, 208896, 274432, 536576, 798720, 1060864, 2109440, 3158016,
  //     81920, 147456, 212992, 278528, 540672, 802816, 1064960, 2113536, 3162112,
  //     98304, 163840, 229376, 294912, 557056, 819200, 1081344, 2129920, 3178496,
  //     114688, 180224, 245760, 311296, 573440, 835584, 1097728, 2146304, 3194880,
  //     327680, 589824, 851968, 1114112, 2162688, 3211264,
  //     393216, 655360, 917504, 1179648, 2228224, 3276800,
  //     458752, 720896, 983040, 1245184, 2293760, 3342336,
  //     1310720, 2359296, 3407872,
  //     1572864, 2621440, 3670016,
  //     1835008, 2883584, 3932160};
  // int mask1[10] = {
  //     1048575, 1048572, 1048560, 1048512, 1048320, 1047552, 1044480, 1032192,
  //     983040, 786432};
  // int mask2[10] = {
  //     3, 15, 63, 255, 1023, 4095, 16383, 65535, 262143, 1048575};


  //
  // Generating mismatch queries and probing the lookup table
  //
  segment_step = word_length + ltable_step - 1;
  //    segment_step = word_length;
  //cout << "Setting segment step to " << segment_step << endl;
  for (shift = 0; shift <= (int)strlen(string) - segment_step; shift = shift + segment_step)
  {
    for (ltable_shift = 0; ltable_shift < ltable_step; ltable_shift++)
    {
      if (shift + ltable_shift + word_length < query_full_len)
      {
        lts = ltable_shift;
      }
      else
      {
        lts = -ltable_shift;
      }
      strncpy(word, &string[shift + lts], word_length);
      //cout << "Searching with " << word << endl;
      for (c = 0; c < MAX_LTABLE; c++)
      {
        visited[c] = 'n';
      }
      //cerr << "SHIFT:" << shift << "+" << lts << endl;
      query_len = strlen(word);
      strncpy(new_word, word,query_len);
      bin = oligo_to_binary(new_word);
      // Generate all the substituted sequences, with 0 to mismatch substitutions
      for (q = 0; q < mismatch_depth[mismatches]; q++)
      {
        if (mismatch_base[q] < Power_of_4[word_length])
        {
          value = bin ^ mismatch_base[q];
          //   cout << "Searching with VALUE=" << value << " MB=" << mismatch_base[q] << " Q=" << q << endl;
          if (visited[value] == 'n')
          {
            for (pos = 1; pos < ltable[value][0]; pos++) //was ltable_count[value]
            {
              //       cout << "End reached. POS=" << pos << " MC =" << match_count << endl;
              matches[match_count] = ltable[value][pos];
              match_shift[match_count] = shift + lts;
              match_count++;
            }
          }
          visited[value] = 'y';
        }
      } // end substituting
      //cout << "Going on to indels..." << endl;
      if (indels)
      {
        // Generate all the deleted sequences, with 1 deletion
        for (q = 2; q < word_length; q++)
        {
          bin1 = bin & mask1[q];
          bin2 = bin & mask2[q - 2];
          bin2 = bin2 << 2;
          for (base = 0; base < 4; base++)
          {
            value = (bin1 | bin2) | base;
            //cout << "Searching with VALUE=" << value << " BIN1=" << bin1 << " MASK1=" << mask1[q] << " BIN2=" << bin2 << " MASK2=" << mask2[q-2] << " Q=" << q << " BASE=" << base << endl;
            // printf("Searching with VALUE=%d BIN1=%d MASK1=%d BIN2=%d MASK2= %d Q= %d BASE=%d\n",value,bin1,mask1[q],bin2,mask2[q-2],q,base); 
            if (visited[value] == 'n')
            {
              for (pos = 1; pos < ltable[value][0]; pos++) //ltable_count
              {
                matches[match_count] = ltable[value][pos];
                match_shift[match_count] = shift + lts;
                match_count++;
              }
            }
          }
          visited[value] = 'y';
        } // end deleting
          // Generate all the inserted sequences, with 1 deletion
        for (q = 2; q < word_length; q++)
        {
          bin1 = bin & mask1[q - 1];
          bin2 = bin & mask2[q - 2];
          bin1 = bin1 << 2;
          bin1 = bin1 & mask2[word_length - 1]; // zero bits shifted too high
          for (base = 0; base < 4; base++)
          {
            value = bin1 | (base << 2 * (q - 1)) | bin2;
            //   cout << "Searching with VALUE=" << value << " BIN1=" << bin1 << " MASK1=" << mask1[q-1] << " BIN2=" << bin2 << " MASK2=" << mask2[q-2] << " Q=" << q << " BASE=" << base << endl;
            if (visited[value] == 'n')
            {
              for (pos = 1; pos < ltable[value][0]; pos++)//ltable_count[value]
              {
                matches[match_count] = ltable[value][pos];
                match_shift[match_count] = shift + lts;
                match_count++;
              }
            }
          }
          visited[value] = 'y';
        } // end inserting

      } // end if(indels)
    }   // end of ltable_shift for loop
  }     // end of shift for loop
}


int oligo_to_binary(char *string)
{
  int j, letter_value, my_value;
  // int Power_of_4[12] = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576};
  // Extract first word_length bases and query the hash
  my_value = 0;
  for (j = 0; j < (int)strlen(string); j++)
  {
    switch (string[j])
    {
    case 'A':
      letter_value = 0;
      break;
    case 'C':
      letter_value = 1;
      break;
    case 'G':
      letter_value = 2;
      break;
    case 'T':
      letter_value = 3;
      break;
    }
    my_value = my_value + letter_value * Power_of_4[j];
  }
  if (my_value > Power_of_4[WORD_LENGTH] - 1)
  {
    //cerr << "Warning: oligo_to_binary returned " << my_value << endl;
    //cerr << "POS: " << global_pos << " OLIGO " << oligo << endl;
  }
  return my_value;
}


void reverse(char* string,char* new_string,size_t length)
{
  char new_base[1];
  // char out[MAX_CHAR];
  strncpy(new_string, "",1);
  for (int i = strlen(string); i > 0; i = i - 1)
  {
    switch (string[i - 1])
    {
    case 'A':
      new_base[0] = 'T';
      break;
    case 'T':
      new_base[0] = 'A';
      break;
    case 'C':
      new_base[0] = 'G';
      break;
    case 'G':
      new_base[0] = 'C';
      break;
    }
    // new_base[1] = '\0';
    strncat(new_string, new_base,1);
  }
  // int len = strlen(new_string);
  char *nl = "\n";
  strncat(new_string,nl,1);
  // new_string = out;
  //return string;
}
