
#include <stdio.h>
#include <stdlib.h>

// provided prototypes
int FillStencil1D(double* stencil, int n_stencil);
int FillStencil2D(double* stencil, int n_stencil);

int FillStencil1D(double* stencil, int n_stencil){

// fill the stencil array
    switch(n_stencil){
        case 3:
            stencil[0] =  1.0;
            stencil[1] = -2.0;
            stencil[2] =  1.0;
            return 0;

        case 5:
            stencil[0] =   -1.0/12.0;
            stencil[1] =   16.0/12.0;
            stencil[2] =  -30.0/12.0;
            stencil[3] =   16.0/12.0;
            stencil[4] =   -1.0/12.0;
            return 0;

        case 7:
            stencil[0] =    2.0/180.0;
            stencil[1] =  -27.0/180.0;
            stencil[2] =  270.0/180.0;
            stencil[3] = -490.0/180.0;
            stencil[4] =  270.0/180.0;
            stencil[5] =  -27.0/180.0;
            stencil[6] =    2.0/180.0;
            return 0;

        case 9:
            stencil[0] =     -9.0/5040.0;
            stencil[1] =    128.0/5040.0;
            stencil[2] =  -1008.0/5040.0;
            stencil[3] =   8064.0/5040.0;
            stencil[4] = -14350.0/5040.0;
            stencil[5] =   8064.0/5040.0;
            stencil[6] =  -1008.0/5040.0;
            stencil[7] =    128.0/5040.0;
            stencil[8] =     -9.0/5040.0;
            return 0;

        case 11:
            stencil[0]  =      8.0/25200.0;
            stencil[1]  =   -125.0/25200.0;
            stencil[2]  =   1000.0/25200.0;
            stencil[3]  =  -6000.0/25200.0;
            stencil[4]  =  42000.0/25200.0;
            stencil[5]  = -73766.0/25200.0;
            stencil[6]  =  42000.0/25200.0;
            stencil[7]  =  -6000.0/25200.0;
            stencil[8]  =   1000.0/25200.0;
            stencil[9]  =   -125.0/25200.0;
            stencil[10] =      8.0/25200.0;
            return 0;

        case 13:
            stencil[0]  =      -50.0/831600.0;
            stencil[1]  =      864.0/831600.0;
            stencil[2]  =    -7425.0/831600.0;
            stencil[3]  =    44000.0/831600.0;
            stencil[4]  =  -222750.0/831600.0;
            stencil[5]  =  1425600.0/831600.0;
            stencil[6]  = -2480478.0/831600.0;
            stencil[7]  =  1425600.0/831600.0;
            stencil[8]  =  -222750.0/831600.0;
            stencil[9]  =    44000.0/831600.0;
            stencil[10] =    -7425.0/831600.0;
            stencil[11] =      864.0/831600.0;
            stencil[12] =      -50.0/831600.0;
            return 0;

        default:
        // if stencil is not implemented print an error message and exit program
            fprintf(stderr,
                "\n (-) Error no data for %d-point stencil available."
                "\n     Aborting - please check your input..."
                "\n\n"
                , n_stencil
            );
            exit(EXIT_FAILURE);
    }
}


int FillStencil2D(double* stencil, int n_stencil){

// fill the stencil array
    switch(n_stencil){
        case 3:
            stencil[0]   =  1.0;
            stencil[1]   =  0.0;
            stencil[2]   =  1.0;

            stencil[3]   =  0.0;
            stencil[4]   =  -4.0;
            stencil[5]   =  0.0;

            stencil[6]   =  1.0;
            stencil[7]   =  0.0;
            stencil[8]   =  1.0;

            return 0;

        case 5:
            stencil[0]   =   0.0;
            stencil[1]   =   0.0;
            stencil[2]   =  -0.166666666666667;
            stencil[3]   =   0.0;
            stencil[4]   =   0.0;

            stencil[5]   =   0.0;
            stencil[6]   =   0.0;
            stencil[7]   =   2.666666666666667;
            stencil[8]   =   0.0;
            stencil[9]   =   0.0;

            stencil[10]  =  -0.166666666666667;
            stencil[11]  =   2.666666666666667;
            stencil[12]  =   -10.0;
            stencil[13]  =   2.666666666666667;
            stencil[14]  =  -0.166666666666667;

            stencil[15]  =   0.0;
            stencil[16]  =   0.0;
            stencil[17]  =   2.666666666666667;
            stencil[18]  =   0.0;
            stencil[19]  =   0.0;

            stencil[20]  =   0.0;
            stencil[21]  =   0.0;
            stencil[22]  =  -0.166666666666667;
            stencil[23]  =   0.0;
            stencil[24]  =   0.0;

            return 0;

        case 7:
            stencil[0]   =    0.0;
            stencil[1]   =    0.0;
            stencil[2]   =    0.0;
            stencil[3]   =    0.022222222222222;
            stencil[4]   =    0.0;
            stencil[5]   =    0.0;
            stencil[6]   =    0.0;

            stencil[7]   =  0.0;
            stencil[8]   =  -0.006944444444444;
            stencil[9]   =  0.027777777777778;
            stencil[10]  =  -0.341666666666667;
            stencil[11]  =  0.027777777777778;
            stencil[12]  =  -0.006944444444444;
            stencil[13]  =  0.0;

            stencil[14]  =  0.0;
            stencil[15]  =  0.027777777777778;
            stencil[16]  =  -0.111111111111111;
            stencil[17]  =  3.166666666666666;
            stencil[18]  =  -0.111111111111111;
            stencil[19]  =  0.027777777777778;
            stencil[20]  =  0.0;

            stencil[21]  =    0.022222222222222;
            stencil[22]  =   -0.341666666666667;
            stencil[23]  =    3.166666666666666;
            stencil[24]  =   -11.138888888888889;
            stencil[25]  =    3.166666666666666;
            stencil[26]  =   -0.341666666666667;
            stencil[27]  =    0.022222222222222;

            stencil[28]  =  0.0;
            stencil[29]  =  0.027777777777778;
            stencil[30]  =  -0.111111111111111;
            stencil[31]  =  3.166666666666666;
            stencil[32]  =  -0.111111111111111;
            stencil[33]  =  0.027777777777778;
            stencil[34]  =  0.0;

            stencil[35]  =  0.0;
            stencil[36]  =  -0.006944444444444;
            stencil[37]  =  0.027777777777778;
            stencil[38]  =  -0.341666666666667;
            stencil[39]  =  0.027777777777778;
            stencil[40]  =  -0.006944444444444;
            stencil[41]  =  0.0;

            stencil[42]  =    0.0;
            stencil[43]  =    0.0;
            stencil[44]  =    0.0;
            stencil[45]  =    0.022222222222222;
            stencil[46]  =    0.0;
            stencil[47]  =    0.0;
            stencil[48]  =    0.0;

            return 0;

        case 9:
            stencil[0]   =   0.0;
            stencil[1]   =   0.0;
            stencil[2]   =   0.0;
            stencil[3]   =   0.0;
            stencil[4]   =  -0.003571428571429;
            stencil[5]   =   0.0;
            stencil[6]   =   0.0;
            stencil[7]   =   0.0;
            stencil[8]   =   0.0;

            stencil[9]   =   0.0;
            stencil[10]  =  -0.000123456790123;
            stencil[11]  =   0.000509259259259;
            stencil[12]  =  -0.000925925925926;
            stencil[13]  =   0.051873897707231;
            stencil[14]  =  -0.000925925925926;
            stencil[15]  =   0.000509259259259;
            stencil[16]  =  -0.000123456790123;
            stencil[17]  =   0.0;

            stencil[18]  =   0.0;
            stencil[19]  =   0.000509259259259;
            stencil[20]  =  -0.001666666666667;
            stencil[21]  =   0.002083333333333;
            stencil[22]  =  -0.401851851851852;
            stencil[23]  =   0.002083333333333;
            stencil[24]  =  -0.001666666666667;
            stencil[25]  =   0.000509259259259;
            stencil[26]  =   0.0;

            stencil[27]  =   0.000000000000000;
            stencil[28]  =  -0.000925925925926;
            stencil[29]  =   0.002083333333333;
            stencil[30]  =   0.0;
            stencil[31]  =   3.197685185185184;
            stencil[32]  =   0.0;
            stencil[33]  =   0.002083333333333;
            stencil[34]  =  -0.000925925925926;
            stencil[35]  =   0.0;

            stencil[36]  =  -0.003571428571429;
            stencil[37]  =   0.051873897707231;
            stencil[38]  =  -0.401851851851852;
            stencil[39]  =   3.197685185185186;
            stencil[40]  =  -11.382716049382715;
            stencil[41]  =   3.197685185185186;
            stencil[42]  =  -0.401851851851852;
            stencil[43]  =   0.051873897707231;
            stencil[44]  =  -0.003571428571429;

            stencil[45]  =   0.000000000000000;
            stencil[46]  =  -0.000925925925926;
            stencil[47]  =   0.002083333333333;
            stencil[48]  =   0.0;
            stencil[49]  =   3.197685185185184;
            stencil[50]  =   0.0;
            stencil[51]  =   0.002083333333333;
            stencil[52]  =  -0.000925925925926;
            stencil[53]  =   0.0;

            stencil[54]  =   0.0;
            stencil[55]  =   0.000509259259259;
            stencil[56]  =  -0.001666666666667;
            stencil[57]  =   0.002083333333333;
            stencil[58]  =  -0.401851851851852;
            stencil[59]  =   0.002083333333333;
            stencil[60]  =  -0.001666666666667;
            stencil[61]  =   0.000509259259259;
            stencil[62]  =   0.0;

            stencil[63]  =   0.0;
            stencil[64]  =  -0.000123456790123;
            stencil[65]  =   0.000509259259259;
            stencil[66]  =  -0.000925925925926;
            stencil[67]  =   0.051873897707231;
            stencil[68]  =  -0.000925925925926;
            stencil[69]  =   0.000509259259259;
            stencil[70]  =  -0.000123456790123;
            stencil[71]  =   0.0;

            stencil[72]  =   0.0;
            stencil[73]  =   0.0;
            stencil[74]  =   0.0;
            stencil[75]  =   0.0;
            stencil[76]  =  -0.003571428571429;
            stencil[77]  =   0.0;
            stencil[78]  =   0.0;
            stencil[79]  =   0.0;
            stencil[80]  =   0.0;

            return 0;

        case 11:
            stencil[0]   =   0.0;
            stencil[1]   =   0.0;
            stencil[2]   =   0.0;
            stencil[3]   =   0.0;
            stencil[4]   =   0.0;
            stencil[5]   =   0.000634920634921;
            stencil[6]   =   0.0;
            stencil[7]   =   0.0;
            stencil[8]   =   0.0;
            stencil[9]   =   0.0;
            stencil[10]  =   0.0;

            stencil[11]  =   0.0;
            stencil[12]  =  -0.000003188775510;
            stencil[13]  =   0.000018345301083;
            stencil[14]  =  -0.000050429894180;
            stencil[15]  =   0.000087632275132;
            stencil[16]  =  -0.010025352733686;
            stencil[17]  =   0.000087632275132;
            stencil[18]  =  -0.000050429894180;
            stencil[19]  =   0.000018345301083;
            stencil[20]  =  -0.000003188775510;
            stencil[21]  =   0.0;

            stencil[22]  =   0.0;
            stencil[23]  =   0.000018345301083;
            stencil[24]  =  -0.000097159234064;
            stencil[25]  =   0.000249118165785;
            stencil[26]  =  -0.000415564373898;
            stencil[27]  =   0.079855599647266;
            stencil[28]  =  -0.000415564373898;
            stencil[29]  =   0.000249118165785;
            stencil[30]  =  -0.000097159234064;
            stencil[31]  =   0.000018345301083;
            stencil[32]  =   0.0;

            stencil[33]  =   0.0;
            stencil[34]  =  -0.000050429894180;
            stencil[35]  =   0.000249118165785;
            stencil[36]  =  -0.000601851851852;
            stencil[37]  =   0.000972222222222;
            stencil[38]  =  -0.477328593474427;
            stencil[39]  =   0.000972222222222;
            stencil[40]  =  -0.000601851851852;
            stencil[41]  =   0.000249118165785;
            stencil[42]  =  -0.000050429894180;
            stencil[43]  =   0.0;

            stencil[44]  =   0.0;
            stencil[45]  =   0.000087632275132;
            stencil[46]  =  -0.000415564373898;
            stencil[47]  =   0.000972222222222;
            stencil[48]  =  -0.001550925925927;
            stencil[49]  =   3.335146604938272;
            stencil[50]  =  -0.001550925925927;
            stencil[51]  =   0.000972222222222;
            stencil[52]  =  -0.000415564373898;
            stencil[53]  =   0.000087632275132;
            stencil[54]  =   0.0;

            stencil[55]  =   0.000634920634921;
            stencil[56]  =  -0.010025352733686;
            stencil[57]  =   0.079855599647266;
            stencil[58]  =  -0.477328593474427;
            stencil[59]  =   3.335146604938272;
            stencil[60]  =  -11.711010802469133;
            stencil[61]  =   3.335146604938272;
            stencil[62]  =  -0.477328593474427;
            stencil[63]  =   0.079855599647266;
            stencil[64]  =  -0.010025352733686;
            stencil[65]  =   0.000634920634921;

            stencil[66]  =   0.0;
            stencil[67]  =   0.000087632275132;
            stencil[68]  =  -0.000415564373898;
            stencil[69]  =   0.000972222222222;
            stencil[70]  =  -0.001550925925927;
            stencil[71]  =   3.335146604938272;
            stencil[72]  =  -0.001550925925927;
            stencil[73]  =   0.000972222222222;
            stencil[74]  =  -0.000415564373898;
            stencil[75]  =   0.000087632275132;
            stencil[76]  =   0.0;

            stencil[77]  =   0.0;
            stencil[78]  =  -0.000050429894180;
            stencil[79]  =   0.000249118165785;
            stencil[80]  =  -0.000601851851852;
            stencil[81]  =   0.000972222222222;
            stencil[82]  =  -0.477328593474427;
            stencil[83]  =   0.000972222222222;
            stencil[84]  =  -0.000601851851852;
            stencil[85]  =   0.000249118165785;
            stencil[86]  =  -0.000050429894180;
            stencil[87]  =   0.0;

            stencil[88]  =   0.0;
            stencil[89]  =   0.000018345301083;
            stencil[90]  =  -0.000097159234064;
            stencil[91]  =   0.000249118165785;
            stencil[92]  =  -0.000415564373898;
            stencil[93]  =   0.079855599647266;
            stencil[94]  =  -0.000415564373898;
            stencil[95]  =   0.000249118165785;
            stencil[96]  =  -0.000097159234064;
            stencil[97]  =   0.000018345301083;
            stencil[98]  =   0.0;

            stencil[99]  =    0.0;
            stencil[100] =  -0.000003188775510;
            stencil[101] =   0.000018345301083;
            stencil[102] =  -0.000050429894180;
            stencil[103] =   0.000087632275132;
            stencil[104] =  -0.010025352733686;
            stencil[105] =   0.000087632275132;
            stencil[106] =  -0.000050429894180;
            stencil[107] =   0.000018345301083;
            stencil[108] =  -0.000003188775510;
            stencil[109] =   0.0;

            stencil[110] =   0.0;
            stencil[111] =   0.0;
            stencil[112] =   0.0;
            stencil[113] =   0.0;
            stencil[114] =   0.0;
            stencil[115] =   0.000634920634921;
            stencil[116] =   0.0;
            stencil[117] =   0.0;
            stencil[118] =   0.0;
            stencil[119] =   0.0;
            stencil[120] =   0.0;

            return 0;

        case 13:
            stencil[0]   =    0.0;
            stencil[1]   =    0.0;
            stencil[2]   =    0.0;
            stencil[3]   =    0.0;
            stencil[4]   =    0.0;
            stencil[5]   =    0.0;
            stencil[6]   =   -0.000120250120250;
            stencil[7]   =    0.0;
            stencil[8]   =    0.0;
            stencil[9]   =    0.0;
            stencil[10]  =    0.0;
            stencil[11]  =    0.0;
            stencil[12]  =    0.0;

            stencil[13]  =    0.0;
            stencil[14]  =   -0.000000100781053;
            stencil[15]  =    0.000000767668178;
            stencil[16]  =   -0.000002791600179;
            stencil[17]  =    0.000006389361300;
            stencil[18]  =   -0.000010196208113;
            stencil[19]  =    0.002094547102417;
            stencil[20]  =   -0.000010196208113;
            stencil[21]  =    0.000006389361300;
            stencil[22]  =   -0.000002791600179;
            stencil[23]  =    0.000000767668178;
            stencil[24]  =   -0.000000100781053;
            stencil[25]  =    0.0;

            stencil[26]  =    0.0;
            stencil[27]  =    0.000000767668178;
            stencil[28]  =   -0.000005691078515;
            stencil[29]  =    0.000020298371406;
            stencil[30]  =   -0.000045892778408;
            stencil[31]  =    0.000072751322751;
            stencil[32]  =   -0.018016014629875;
            stencil[33]  =    0.000072751322751;
            stencil[34]  =   -0.000045892778408;
            stencil[35]  =    0.000020298371406;
            stencil[36]  =   -0.000005691078515;
            stencil[37]  =    0.000000767668178;
            stencil[38]  =    0.0;

            stencil[39]  =    0.0;
            stencil[40]  =   -0.000002791600179;
            stencil[41]  =    0.000020298371406;
            stencil[42]  =   -0.000071570294785;
            stencil[43]  =    0.000160967550600;
            stencil[44]  =   -0.000254721487360;
            stencil[45]  =    0.106710978835979;
            stencil[46]  =   -0.000254721487360;
            stencil[47]  =    0.000160967550600;
            stencil[48]  =   -0.000071570294785;
            stencil[49]  =    0.000020298371406;
            stencil[50]  =   -0.000002791600179;
            stencil[51]  =    0.000000000000000;

            stencil[52]  =    0.0;
            stencil[53]  =    0.000006389361300;
            stencil[54]  =   -0.000045892778408;
            stencil[55]  =    0.000160967550600;
            stencil[56]  =   -0.000361906336609;
            stencil[57]  =    0.000573467813051;
            stencil[58]  =   -0.539951765505585;
            stencil[59]  =    0.000573467813051;
            stencil[60]  =   -0.000361906336609;
            stencil[61]  =    0.000160967550600;
            stencil[62]  =   -0.000045892778408;
            stencil[63]  =    0.000006389361300;
            stencil[64]  =    0.0;

            stencil[65]  =    0.0;
            stencil[66]  =   -0.000010196208113;
            stencil[67]  =    0.000072751322751;
            stencil[68]  =   -0.000254721487360;
            stencil[69]  =    0.000573467813051;
            stencil[70]  =   -0.000910493827160;
            stencil[71]  =    3.454629813345090;
            stencil[72]  =   -0.000910493827160;
            stencil[73]  =    0.000573467813051;
            stencil[74]  =   -0.000254721487360;
            stencil[75]  =    0.000072751322751;
            stencil[76]  =   -0.000010196208113;
            stencil[77]  =    0.0;

            stencil[78]  =   -0.000120250120250;
            stencil[79]  =    0.002094547102417;
            stencil[80]  =   -0.018016014629875;
            stencil[81]  =    0.106710978835979;
            stencil[82]  =   -0.539951765505585;
            stencil[83]  =    3.454629813345090;
            stencil[84]  =   -12.020383506944441;
            stencil[85]  =    3.454629813345090;
            stencil[86]  =   -0.539951765505585;
            stencil[87]  =    0.106710978835979;
            stencil[88]  =   -0.018016014629875;
            stencil[89]  =    0.002094547102417;
            stencil[90]  =   -0.000120250120250;

            stencil[91]  =    0.0;
            stencil[92]  =   -0.000010196208113;
            stencil[93]  =    0.000072751322751;
            stencil[94]  =   -0.000254721487360;
            stencil[95]  =    0.000573467813051;
            stencil[96]  =   -0.000910493827160;
            stencil[97]  =    3.454629813345090;
            stencil[98]  =   -0.000910493827160;
            stencil[99]  =    0.000573467813051;
            stencil[100] =   -0.000254721487360;
            stencil[101] =    0.000072751322751;
            stencil[102] =   -0.000010196208113;
            stencil[103] =    0.0;

            stencil[104] =    0.0;
            stencil[105] =    0.000006389361300;
            stencil[106] =   -0.000045892778408;
            stencil[107] =    0.000160967550600;
            stencil[108] =   -0.000361906336609;
            stencil[109] =    0.000573467813051;
            stencil[110] =   -0.539951765505585;
            stencil[111] =    0.000573467813051;
            stencil[112] =   -0.000361906336609;
            stencil[113] =    0.000160967550600;
            stencil[114] =   -0.000045892778408;
            stencil[115] =    0.000006389361300;
            stencil[116] =    0.0;

            stencil[117] =    0.0;
            stencil[118] =   -0.000002791600179;
            stencil[119] =    0.000020298371406;
            stencil[120] =   -0.000071570294785;
            stencil[121] =    0.000160967550600;
            stencil[122] =   -0.000254721487360;
            stencil[123] =    0.106710978835979;
            stencil[124] =   -0.000254721487360;
            stencil[125] =    0.000160967550600;
            stencil[126] =   -0.000071570294785;
            stencil[127] =    0.000020298371406;
            stencil[128] =   -0.000002791600179;
            stencil[129] =    0.000000000000000;

            stencil[130] =    0.0;
            stencil[131] =    0.000000767668178;
            stencil[132] =   -0.000005691078515;
            stencil[133] =    0.000020298371406;
            stencil[134] =   -0.000045892778408;
            stencil[135] =    0.000072751322751;
            stencil[136] =   -0.018016014629875;
            stencil[137] =    0.000072751322751;
            stencil[138] =   -0.000045892778408;
            stencil[139] =    0.000020298371406;
            stencil[140] =   -0.000005691078515;
            stencil[141] =    0.000000767668178;
            stencil[142] =    0.0;

            stencil[143] =    0.0;
            stencil[144] =   -0.000000100781053;
            stencil[145] =    0.000000767668178;
            stencil[146] =   -0.000002791600179;
            stencil[147] =    0.000006389361300;
            stencil[148] =   -0.000010196208113;
            stencil[149] =    0.002094547102417;
            stencil[150] =   -0.000010196208113;
            stencil[151] =    0.000006389361300;
            stencil[152] =   -0.000002791600179;
            stencil[153] =    0.000000767668178;
            stencil[154] =   -0.000000100781053;
            stencil[155] =    0.0;

            stencil[156] =    0.0;
            stencil[157] =    0.0;
            stencil[158] =    0.0;
            stencil[159] =    0.0;
            stencil[160] =    0.0;
            stencil[161] =    0.0;
            stencil[162] =   -0.000120250120250;
            stencil[163] =    0.0;
            stencil[164] =    0.0;
            stencil[165] =    0.0;
            stencil[166] =    0.0;
            stencil[167] =    0.0;
            stencil[168] =    0.0;

            return 0;

        default:
        // if stencil is not implemented print an error message and exit program
            fprintf(stderr,
                "\n (-) Error no data for %d-point stencil available."
                "\n     Aborting - please check your input..."
                "\n\n"
                , n_stencil
            );
            exit(EXIT_FAILURE);
    }
}
