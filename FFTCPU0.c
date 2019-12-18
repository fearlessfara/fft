/*
 ============================================================================
 Name        : FFTCPU0.c
 Author      : Faraone Christian Gennaro & Professor Michael Ryan
 Version     : 1.0
 Copyright   : Copyright OCE Technology
 Description : Multicore FFT

 This program demonstrates the use of 1, 2, or 4 CPUs on a S698PM
 multi-core processor to calculate a Fast Fourier Transform (FFT).

 A selection of different input arrays is provided,
 additional input arrays can be created easily.

 At present the input array length must be a power of 2,
 and the scalar type must be C float.

 The functions used are defined in fft.h and documented there.
 Please consult fft.h for more information

 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include "fft.h"

/*
 * main()
 *
 * Usage:
 *
 * 		Assign values to each of the variables below
 *
 * 			lenP2	 		the power of 2 that gives the length of the input array
 *							(at present from 5 to 17, i.e. length 32 to 131,072)
 *
 * 			whichcpus		the CPUs to use, from 1 to 0xf (binary xxxx, 1 = use)
 *
 * 			numruns			number of test runs - timings are given for total runs
 *
 * 			type 			0 for forward FFT, 1 for inverse FFT
 *
 * 			normalise		0 for not normalized, 1 for normalized
 *
 * 			arrange			0 if data already arranged, 1 if need to arrange
 *
 *    		realonly		0 for (in-phase, quadrature) input, 1 for (in-phase, all 0);
 *
 * 			displayLines	number of lines of data to display
 *
 * 		Set up the complex data input in cin.
 * 		(a number of choices are already set up and commented out, it is simple to add more)
 *
 *		Compile the program.
 *
 *		Use the DMON script provided to load and run the program
 *		This also sets up the CPU stacks and enables the L2 cache.
 *
 *		N.B. L2 cache is disabled by default,
 *			 if not switched on timings are affected dramatically.
 *
 */
int main(int argc,char ** argv)
{
	/*
	 * Initial entry point for all CPUs  -- DO NOT CHANGE
	 */
	if(START_CPU_INDEX != fftCpusStart())		// N.B. only START_CPU returns from this call
	{
		printf("\nAbandoning program, startup problem\n");
		return 0;
	}	// if

	/*
	 * START HERE	- select input size, CPUs to use, number of test runs
	 */
	int lenP2 = 12;					// input array length is this power of 2 (from 5 to 17 at present)

	int whichcpus = 0xf;			// CPUs to use, hex value xxxx, x = 1 to turn on, 0xf for all on

    const int numruns = 200;		// number of test runs

    const int type = 0;				// 0 for forward FFT, 1 for inverse FFT

    const int normalise = 0;		// 0 for not normalized, 1 for normalized

    const int arrange = 1;			// 0 for don't arrange, 1 for arrange

    const int realonly = 0;			// 0 for (in-phase, quadrature), 1 for (in-phase,all 0);

    int displayLines = 64;			// number of lines of the data to display, if enabled below

    // check input and correct if possible, give message and exit if settings incorrect
    if(settingsIncorrect(lenP2, &whichcpus, numruns, type, &displayLines)) return 0;

	/*
	 *  VARIOUS EXAMPLES
	 *  Example lengths depend on lenP2 above
	 *  The examples not in use are commented out
	 *  Other input choices can be set up as (in-phase, quadrature)in cin as desired
	 *  N.B. If realonly is not selected both components should be set up,
	 *  	 otherwise previous data may cause anomalous results.
	 *  	 If realonly is selected, the quadrature part of the input is overwritten
	 *  	 with 0s.
	 */
	fft_cpx * cin = (void *)FFT_CIN;	// start of complex number input data array
	int length = 1 << lenP2;			// length of complex numbers input array

    int i;

	// test: combination of frequencies and phases
	//for (i=0;i<length;++i) {
//		cin[i].re = 1.0 + sin(M_PI*i/(length >> 1))			// in-phase part
	//				    + 3*cos(M_PI*2*i/(length >> 1))
	//				    + 2*cos(M_PI*6*i/(length >> 1))
	//					+ 7*cos(M_PI*((length >> 1) - 1)*i/(length >> 1))
	//				    + 5*cos(M_PI*i);					// Nyquist
	//    cin[i].im = sin(M_PI*3*i/(length >> 1));			// quadrature part
//	}	// for

	// test: ramp
//	for (i=0;i<length;++i) {
//		cin[i].re = i;
//	    cin[i].im = 0;
//	}	// for

	// test: random
//	srand(time(0));
//	for (i=0;i<length;++i) {
//      cin[i].re = rand();
//      cin[i].im = 0;
//  }	// for
// square wave
	for (i=0;i<(length >> 1);++i) {
      cin[i].re = 3;
      cin[i].im = 0;
  }	// for

  	for (i=(length >> 1); i < length;++i) {
      cin[i].re = 0;
      cin[i].im = 0;
  }	// for
  
  
	// show the input data
//	printf("Complex data:\n");
//	showC(cin, displayLines, 0);

	/*
	 * Do FFT test
	 * Does numruns, first with one processor, then with the requested processors (whichcpus)
	 */
	// set up the twiddles array - its length is half that of the input array.
	fft_cpx * tw = (void *)FFT_TWIDDLES;	// preassigned space in SRAM starting at this address
	fillTwiddles(tw,lenP2,type);			// 0 for forward FFT, 1 for inverse FFT

	// if data is in-phase only, i.e. no quadrature component, pack into first half
	if(realonly)
	{
		fftRealOnlyArrange(cin,lenP2);
		lenP2 -= 1;							// half the size, but leave value of length as is
	}	// if

	// if need to arrange, set up mapping array - N.B. it is its own inverse
	if(arrange)
	{
		fftArrangeIndex((void *)FFT_MAP, lenP2);
	}	// if

	// set up a pointer for rearranged input and clear it
	// if 'arrange' set result will be here
	fft_cpx * rain = (void *)FFT_RAIN; 		// preassigned space in SRAM starting at this address
	for(i = 0; i < length; i++)
	{
		rain[i].re = 0;
		rain[i].im = 0;
	}

//	int map[1 << lenP2];
//	fftArrangeIndex(map, lenP2);
//	for (i=0;i<length;++i) {
//		rain[map[i]] = cin[i];
//	}	// for

	// start the cpus being used other than this one - they will halt waiting for work
	startCpus(whichcpus);

	// get ready to time
	clock_t start, finish;					// start and finish times
	unsigned long single, multiple;			// times with single, multiple C{Us

	printf("\n\nDoing %d runs of %d length FFT with 1 processor\n", numruns, length);

	start = clock();
	for(i = 0; i < numruns; i++)
	{
		if(arrange)							// will do rearrangement in parallel
		{
			// do it using this CPU, the startup CPU, usually CPU 0
			fftMulti(cin, lenP2, tw, 1<<START_CPU_INDEX, normalise, arrange, realonly);
		} else {

			// first, rearrange the input into rain
			fftArrangeOut(cin, rain,lenP2);

			// do it using this CPU, the startup CPU, usually CPU 0
			fftMulti(rain, lenP2, tw, 1<<START_CPU_INDEX, normalise, arrange, realonly);
		}	// else
	}	// for

	finish = clock();
	single = finish - start;

	printf("\n\nDone %d runs of %d length FFT with 1 processor, time: %ld\n", numruns, length, single);

	int numcpus = getNumCpus(whichcpus);	// how many CPUs being used

	printf("\nDoing %d runs of %d length FFT with %d processors\n", numruns, length, numcpus);

	start = clock();
	for(i = 0; i < numruns; i++)
	{

		if(arrange)						// will do rearrangement in parallel
		{
			// do it using selected cpus
			fftMulti(cin, lenP2, tw, whichcpus, normalise, arrange, realonly);
		} else {
			// first, rearrange the input
			fftArrangeOut(cin, rain, lenP2);
			fftMulti(rain, lenP2, tw, whichcpus, normalise, arrange, realonly);
		}	// else
	}	// for

	finish = clock();
	multiple = finish - start;

	/*
	 * FINISH - show FFT result and times
	 */
	printf("Complex FFT:\n");
	show8C(rain, length, 0.001, realonly);

	printf("\n\n%d runs of %d length FFT with %d processor, time: %ld\n", numruns, length, 1, single);
	printf("\n%d runs of %d length FFT with %d processors, time: %ld\n", numruns, length, numcpus, multiple);
	if(realonly)
	{
		printf("\nSpeedup factor for %d real data with %d processors: %2.2f\n", length, numcpus, ((float)single)/(float)multiple);
	} else {
		printf("\nSpeedup factor for %d (in-phase,quadrature) with %d processors: %2.2f\n", length, numcpus, ((float)single)/(float)multiple);
	}	// else

	return 0;

}	// main()


// see fft.h for description
static inline unsigned int get_asr17(void)
{
  unsigned int reg;
  __asm__ (" mov %%asr17, %0 " : "=r"(reg) :);
  return reg;
}


// see fft.h for description
int fftCpusStart(void)
{
	volatile unsigned int status = 0;
	volatile fft_cpu *thisCpu;

	// find out which CPU is running
	int LEON3_Cpu_Index = (get_asr17() >> 28) & 3;		// max 4 cpus, 0 .. 3

	// set up reference to data for this CPU
	thisCpu = (void *)(FFT_FLAGS + LEON3_Cpu_Index*sizeof(fft_cpu));

	// if start CPU (usually CPU0), initialize flags area
	if (START_CPU_INDEX == LEON3_Cpu_Index){
		// START_CPU (usually 0) initializes settings and returns
		volatile fft_cpu *thisCpu;
		int i = 0;
		for(i = 0; i < MAX_CPUS; i++)
		{
			thisCpu = (void *)(FFT_FLAGS + i*sizeof(fft_cpu));

			thisCpu->status = 0;
			thisCpu->cinadrs = 0;
			thisCpu->lenP2 = 0;
			thisCpu->begin = 0;
			thisCpu->nextbegin = 0;
			thisCpu->oddadd = 0;
			thisCpu->startP2 = 0;
			thisCpu->twadrs = 0;
			thisCpu->toffset = 0;
			thisCpu->tinc = 0;
			thisCpu->normalise = 0;
			thisCpu->arrange = 0;
			thisCpu->realonly = 0;
			thisCpu->spare1 = 0x01010101;
			thisCpu->spare2 = 0x10101010;
		}	// for

	}else {
		// other CPUS loop forever, check for work, do it, set done status
		while(1)
		{
			thisCpu->spare1 = 0xcafe;						// for debug

			status = thisCpu->status;
			if(FFT_START == (FFT_WORK_MASK & status))
			{
				thisCpu->spare1 = 0xbabecafe;				// for debug

				fftProc((void *)thisCpu->cinadrs, thisCpu->lenP2,
						thisCpu->begin, thisCpu->nextbegin,
						thisCpu->oddadd, thisCpu->startP2,
						(void *)thisCpu->twadrs,
						thisCpu->toffset, thisCpu->tinc,
						thisCpu->normalise,thisCpu->arrange,
						thisCpu->realonly, thisCpu->secondpass);

				thisCpu->spare2 = 0xcafebabe;				// for debug
				thisCpu->status &= ~FFT_WORK_MASK;			// finished

			}	// if

			// halt this CPU, will be restarted when work available
			__asm__ (" wr %g0, %asr19 ");

		}	// while

	}	// else

	return LEON3_Cpu_Index;		// should only get here for CPU 0

}	// fftCpusStart()


// see fft.h for description
void startCpus(const int whichcpus)
{
	int i;
	int LEON3_Cpu_Index = (get_asr17() >> 28) & 3;		// max 4 cpus, 0 .. 3
	int this_cpu = 1 << LEON3_Cpu_Index;

	int current_cpu = LAST_CPU;							// usually 8

	while(current_cpu)
	{
		for(i=0; i < CPU_START_PAUSE; i++){}						// pause
		if((current_cpu & whichcpus)&&(current_cpu != this_cpu))
		{
			*((int *)IRQMP_ADRS) = 0x380b0000 | current_cpu;
		}	// if

		current_cpu >>= 1;

	}	// while

	for(i=0; i < CPU_START_PAUSE; i++){}							// pause
	return;
}	// startCpus()


/*
 * INPUT ARRAY REARRANGEMENT
 */

// see fft.h for description
void fftArrangeOut(fft_cpx* cin, fft_cpx* result, int lenP2)
{

	// validation check
	if((FFT_MIN_LENP2 > lenP2)||(FFT_MAX_LENP2 < lenP2))
	{
		printf("\nInvalid power of 2 for array length: %d, should be between %d and %d\n",
				lenP2, FFT_MIN_LENP2, FFT_MAX_LENP2);
		return;
	}	// if

	int len = 1 << lenP2;
	int mask = len -1;

	int i = 0;
    for(i = 0; i < len; i++){

    	/*
    	 * Find the index to swap with.
    	 * - treating index as potentially 32 bits,
    	 *   number of bits in use should be less
    	 */
    	int v = i;
    	// swap odd and even bits
    	v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
    	// swap consecutive pairs
    	v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
    	// swap nibbles ...
    	v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
    	// swap bytes
    	v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
    	// swap 2-byte long pairs
    	v = ( v >> 16             ) | ( v               << 16);

    	v >>= (32 - lenP2);				// N.B. arithmetic shift
    	v &= mask;

	    if(v >= i)					// otherwise already done
	    {
	    	result[i] = cin[v];
	    	result[v] = cin[i];
	    }	// if

    }	// for

    return;

}	// fftArrangeOut()


// see fft.h for description
void fftArrangeIn(fft_cpx* x, int lenP2)
{
	// validation check
	if((FFT_MIN_LENP2 > lenP2)||(FFT_MAX_LENP2 < lenP2))
	{
		printf("\nInvalid power of 2 for array length: %d, should be between %d and %d\n",
				lenP2, FFT_MIN_LENP2, FFT_MAX_LENP2);
		return;
	}	// if

	// do it
	int len = 1 << lenP2;			// power of 2
	int mask = len - 1;

	fft_cpx temp;					// for swap

	int i = 0;
    for(i = 0; i < len; i++){		// must check full array, not just half

    	// Find the index to swap with
    	// - treating index as no more than 32 bits,
    	int v = i;
    	// swap odd and even bits
    	v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
    	// swap consecutive pairs
    	v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
    	// swap nibbles ...
    	v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
    	// swap bytes
    	v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
    	// swap 2-byte long pairs
    	v = ( v >> 16             ) | ( v               << 16);

    	v >>= (32 - lenP2);			// arithmetic shift
    	v &=  mask;					// kill off any leading 1s due to arithmetic shift

	    if(v > i)					// otherwise already done or same
	    {
	    	temp = x[i];
	    	x[i] = x[v];
	    	x[v] = temp;
	    }	// if

    }	// for

    return;

}	// fftArrangeIn()


// see fft.h for description
void fftArrangeIndex(int* mapindex, int lenP2)
{
	// validation check
	if((FFT_MIN_LENP2 > lenP2)||(FFT_MAX_LENP2 < lenP2))
	{
		printf("\nInvalid power of 2 for array length: %d, should be between %d and %d\n",
				lenP2, FFT_MIN_LENP2, FFT_MAX_LENP2);
		return;
	}	// if

	// do it
	int len = 1 << lenP2;			// power of 2
	int mask = len - 1;

	int i = 0;
    for(i = 0; i < len; i++){		// must check full array, not just half

    	// Find the index to swap with
    	// - treating index as no more than 32 bits,
    	int v = i;
    	// swap odd and even bits
    	v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
    	// swap consecutive pairs
    	v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
    	// swap nibbles ...
    	v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
    	// swap bytes
    	v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
    	// swap 2-byte long pairs
    	v = ( v >> 16             ) | ( v               << 16);

    	v >>= (32 - lenP2);			// arithmetic shift
    	v &=  mask;						// kill off any leading ones

    	mapindex[i] = v;

    }	// for

    return;

}	// fftArrangeIndex()


// see fft.h for description
void fftRealOnlyArrange(fft_cpx * cin, int lenP2)
{
	int i,j;
	int halflen = 1 << (lenP2-1);

	for (i = 0; i < halflen; i++)
	{
		j = i << 1;
		cin[i].re = cin[j].re;
		cin[i].im = cin[j+1].re;
	}	// for

	for (i=halflen; i < (1 << lenP2); i++)
	{
		cin[i].re = 0;
		cin[i].im = 0;
	}	// for

}	// fftRealOnlyArrange()


/*
 * TWIDDLES ARRAY
 */

// see fft.h for description
void fillTwiddles(fft_cpx* twiddles, const int lenP2, const int type)
{
	// validation checks
	if((FFT_MIN_LENP2 > lenP2)||(FFT_MAX_LENP2 < lenP2))
	{
		printf("\nfillTwiddles(): Invalid power of 2 for array length: %d, should be between %d and %d\n",
				lenP2, FFT_MIN_LENP2, FFT_MAX_LENP2);
		return;
	}	// if

    if ((0 != type)&& (1 != type))		// invalid, should be forward(0) or inverse(1)
    {
		printf("\nInvalid type for twiddles array, should be 0 for forward, 1 for inverse, actual was: 0x%x\n",
				type);
		return;
    }	// if

    /*
     * Do it - using doubles
     */
    int n = 1 << lenP2;
    double twoPiDivNsigned;

    if(1 != type)						// negative exponent for forward transform, positive for inverse
    {
        twoPiDivNsigned= -M_PI / (double)(n >> 1);  // 2*Pi/n
    } else {
        twoPiDivNsigned = M_PI / (double)(n >> 1);	// positive exponent for inverse transform
    }	// else

    double root2Inv  = 1.0/sqrt(2.0);	// quicker than sine or cosine

    /*
     * Fill in the values
     */
    int quarter;						// will be N/4, quarter way round
    int eight;							// will be N/8, eight of way around

    twiddles[0].re = 1.0;				// always same
    twiddles[0].im = 0;

    if(1 != type)						// forward FFT - N.B. using negative exponent
    {
    	quarter = n >> 2;
        twiddles[quarter].re = 0;		// quarter way around with negative exponent (0,-i)
        twiddles[quarter].im = -1;

        eight = quarter >> 1;			// eight of way around (1/root2, -i*1/root2)
        twiddles[eight].re = root2Inv;
        twiddles[eight].im = -root2Inv;

        twiddles[quarter+eight].re = -root2Inv;	// three eights of way around (-1/root2, -i*1/root2)
        twiddles[quarter+eight].im = -root2Inv;

        int j;
        for(j=1; j<eight; j++)
        {
        	double theta = (double)j * twoPiDivNsigned;	// negative

        	twiddles[j].re = cos(theta);			// positive
        	twiddles[j].im = sin(theta);			// negative

        	int pos = quarter-j;
        	twiddles[pos].re = -twiddles[j].im;
        	twiddles[pos].im = -twiddles[j].re;

        	pos = quarter+j;
           	twiddles[pos].re =  twiddles[j].im;
            twiddles[pos].im = -twiddles[j].re;

        	pos = (quarter << 1) - j;
        	twiddles[pos].re = -twiddles[j].re;
        	twiddles[pos].im =  twiddles[j].im;
        }	// for

    } else {										// inverse FFT - N.B. uses positive exponent

        quarter = n >> 2;
        twiddles[quarter].re = 0;					// quarter way around (0,i)
        twiddles[quarter].im = 1;					// N.B. positive exponent

        eight = quarter >> 1;						// eight of way around (1/root2, i*1/root2)
        twiddles[eight].re = root2Inv;
        twiddles[eight].im = root2Inv;

        twiddles[quarter+eight].re = -root2Inv;		// three eights of way around (-1/root2, i*1/root2)
        twiddles[quarter+eight].im =  root2Inv;

        int j;
        for(j=1; j<eight; j++)
        {
        	double theta = (double)j * twoPiDivNsigned;		// positive

        	twiddles[j].re = cos(theta);			// positive
        	twiddles[j].im = sin(theta);			// positive

        	int pos = quarter-j;
        	twiddles[pos].re = twiddles[j].im;
        	twiddles[pos].im = twiddles[j].re;

        	pos = quarter+j;
        	twiddles[pos].re =  -twiddles[j].im;
        	twiddles[pos].im = twiddles[j].re;

        	pos = (quarter << 1) - j;
        	twiddles[pos].re = -twiddles[j].re;
        	twiddles[pos].im =  twiddles[j].im;
        }	// for

    }	// else

}	// fillTwiddles()


/*
 * FFT PROCESSING
 */

// see fft.h for description
void fftProc(fft_cpx* cin, const int lenP2,
			 const int start, const int nextstart,
			 const int oddadd, const int startP2,
    		 const fft_cpx* twiddles,
    		 const int toffset, const int twinc,
    		 const int normalise, const int arrange,
    		 const int realonly, const int secondpass)
{
	/*
	 * Input validation
	 */
	if((FFT_MIN_LENP2 > lenP2)||(FFT_MAX_LENP2 < lenP2))
	{
		printf("\nfftProc():Invalid power of 2 for array length: %d, should be between %d and %d\n",
				lenP2, FFT_MIN_LENP2, FFT_MAX_LENP2);
		return;
	}	// if

	int length = 1 << lenP2;
	int subarraylen = nextstart - start;

    if((0 > start)||(0 > nextstart)||(0 >= subarraylen))
    {
    	printf("\nInvalid input sub-array indices, start index: %d, next start: %d, sub-array length %d\n", start, nextstart, subarraylen);
    	return;
    }	// if

    if ((start >= length)||(nextstart > length))
    {
    	printf("\nSub-array indices incompatible with array length, start index: %d, next start: %d, array length: %d\n",
    			start, nextstart, length);
    	return;
    }	// if

    if (subarraylen & (subarraylen -1))	// non-zero if not power of 2
    {
    	printf("\nInvalid sub-array length, not a power of 2, start index: %d, next start: %d, sub-array length: %d\n",
    			start, nextstart, subarraylen);
    	return;
    }	// if
//printf("fftProc(): start %d, nextstart %d, oddadd %d, startP2 %d, toffset %d, tinc %d, realonly %d\n",
//			start,nextstart, oddadd, startP2, toffset, twinc, realonly);

	/*
	 * Processing
	 */
    int bi;					// index where this butterfly starts
    int ei;					// index of even
    int oi;					// index of odd

    fft_cpx temp;			// holds current twiddles value
    int tstart = toffset;	// position in twiddles array from which to begin
    int tindex;				// index into twiddles array
    int tlimit;				// tindex should be less than this
    int tinc = twinc;		// increment used stepping through twiddles array
    int lasttinc = 1;		// when tinc gets to this, on last butterfly

    int oddoffset;			// amount to add to even index to get odd index in inputarray
    int butterflysize;		// holds current butterfly size
    int bsh;    			// half of the butterfly size

    float norm;				// used to hold normalizing divisor
    register float twre;	// holds real part of twiddles entry
    register float twim;	// holds imaginary part of twiddles entry

    fft_cpx * arranged = cin;
    fft_cpx * rain = (void *)FFT_RAIN;
    int * map = (void *)FFT_MAP;

    /*
     * If this is the second pass, do post processing of FFT result obtained so far
     * N.B. the FFT data being processed is in the second half of the underlying array,
     * 		the FFT result is put in the first half of the underlying array plus the Nyquist
     * 		component at index [length], the remaining entries in the underlying array
     * 		are not updated with negative frequencies, these are just the complex conjugates
     * 		of the positive frequencies.
     */
    if(secondpass)
    {
    	float value;
    	int i = start;
    	int j = (length << 1) - start;

    	// if at start of twiddles array, avoid need to use value at index [2*length]
    	if(0 == start)
    	{
    		value =  (arranged[0].re + arranged[length].re
    		    	  + arranged[0].im + arranged[length].im)/2.0;

    		arranged[0].im = (arranged[0].im - arranged[length].im
    		     			  - arranged[0].re + arranged[length].re)/2.0;

    		arranged[0].re = value;

    		i++;
    		j--;
    	}

    	// do the remaining entries
    	for(; i < nextstart; i++)
    	{
    		value =  (arranged[i].re + arranged[j].re
    				  + twiddles[i].re*(arranged[i].im + arranged[j].im)
    				  + twiddles[i].im*(arranged[i].re - arranged[j].re))/2.0;

    	    arranged[i].im = (arranged[i].im - arranged[j].im
     			         	  + twiddles[i].im*(arranged[i].im + arranged[j].im)
     			         	  - twiddles[i].re*(arranged[i].re - arranged[j].re))/2.0;

    	    arranged[i].re = value;

    	    j--;
    	}	// for

    	// Nyquist case
       	if(length == nextstart)
    	{
       		arranged[length].re = arranged[length].re - arranged[length].im;
       		arranged[length].im = 0;
    	}	// if

    	return;
    }

    /*
     * Do FFT - only get here if this is not second pass
     *
     * N.B. If the original data was in-phase only,
     * 		the underlying array is assumed to be 2*length
     */

    // if real only is true, adjust the twiddle array steps
    if(realonly)
    {
    	tstart <<= 1;		// twiddles was set up for full-length complex array
    	tinc <<= 1;			// rather than for half length reals only
    	lasttinc <<= 1;
    }	// if

    if(arrange)
    {
    	int i;
    	for(i = start; i < nextstart; i++)
    	{
    		rain[i] = arranged[map[i]];
    	}	// for

    	arranged = rain;
    }	// if

    // double size of butterfly each time around outer loop for this sub-array
    for(butterflysize = 1 << startP2; butterflysize <= subarraylen; butterflysize <<= 1)
    {
     	bsh = butterflysize >> 1;
     	oddoffset = bsh + oddadd;

     	if(lasttinc != tinc)	// not last butterfly
     	{
     		// do all the butterflies in this sub-array for this size of butterfly
     		for( bi = start; bi < nextstart; bi += butterflysize)
     		{
     			ei = bi;
     			oi = ei + oddoffset;

     			tindex = tstart;
     			tlimit = tstart + bsh*tinc;

     			// as log n repeats, worth doing first case separately
     	    	if(0 == tindex)
     	    	{
     				tindex += tinc;

     	    		temp.re = arranged[oi].re;
     				temp.im = arranged[oi].im;

     				arranged[oi].re = arranged[ei].re - temp.re;
     				arranged[oi++].im = arranged[ei].im - temp.im;
     				arranged[ei].re += temp.re;
     				arranged[ei++].im += temp.im;
     	    	}	// if

     			// doing each butterfly
     			while(tindex < tlimit)
     			{
     				twre = twiddles[tindex].re;
     				twim = twiddles[tindex].im;

     				tindex += tinc;

     				temp.re = arranged[oi].re*twre
     						  - arranged[oi].im*twim;
     				temp.im = arranged[oi].re*twim
     				     	  + arranged[oi].im*twre;

     				arranged[oi].re = arranged[ei].re - temp.re;
     				arranged[oi++].im = arranged[ei].im - temp.im;
     				arranged[ei].re += temp.re;
     				arranged[ei++].im += temp.im;




     			}	// while

     		}	// for

     	} else {	// last butterfly, only one butterfly to do. Normalize if requested.

     		ei = start;
     		oi = ei + oddoffset;

     		tindex = tstart;
     		tlimit = tindex + lasttinc*bsh;		// lasttinc will be 1 or 2

     		if(normalise)
     		{
     			norm = (float)(realonly ?(1 << (lenP2 + 1)): 1 << lenP2);

     			// not worth doing single case separately
//    	    	if(0 == tindex)
//     	    	{
//     				tindex += tinc;
//
//     	    		temp.re = arranged[oi].re;
//     				temp.im = arranged[oi].im;
//
//     				arranged[oi].re = (arranged[ei].re - temp.re)/norm;
//     				arranged[oi++].im = (arranged[ei].im - temp.im)/norm;
//     				arranged[ei].re += temp.re;
//     				arranged[ei].re /= 2.0;
//     				arranged[ei].im += temp.im;
//     				arranged[ei++].im /= 2.0;
//     	    	}	// if

     			while(tindex < tlimit)
     			{
//     				twre = twiddles[tindex].re;		// not faster
//     				twim = twiddles[tindex].im;
//
//     				tindex += lasttinc;
//
//     				temp.re = arranged[oi].re*twre
//     						  - arranged[oi].im*twim;
//     				temp.im = arranged[oi].re*twim
//     				     	  + arranged[oi].im*twre;

     				temp.re =  arranged[oi].re*twiddles[tindex].re
     						   - arranged[oi].im*twiddles[tindex].im;
     				temp.im =  arranged[oi].re*twiddles[tindex].im
     						   + arranged[oi].im*twiddles[tindex].re;

     				tindex += lasttinc;

     				arranged[oi].re = (arranged[ei].re - temp.re)/norm;
     				arranged[oi++].im = (arranged[ei].im - temp.im)/norm;
     				arranged[ei].re += temp.re;
     				arranged[ei].re /= norm;
     				arranged[ei].im += temp.im;
     				arranged[ei++].im /= norm;

     			}	// while

     		} else {		// not normalizing

     			// not worth doing single case separately
//     			if(0 == tindex)
//     			{
//     				tindex += tinc;
//
//     		   		temp.re = arranged[oi].re;
//     			    temp.im = arranged[oi].im;
//
//     			    arranged[oi].re = arranged[ei].re - temp.re;
//     			    arranged[oi++].im = arranged[ei].im - temp.im;
//     			    arranged[ei].re += temp.re;
//     			    arranged[ei++].im += temp.im;
//     			}	// if

     			while(tindex < tlimit)
     			{
//     				twre = twiddles[tindex].re;	// not faster
//     				twim = twiddles[tindex].im;
//
//     				tindex += lasttinc;
//
//     				temp.re = arranged[oi].re*twre
//     						  - arranged[oi].im*twim;
//     				temp.im = arranged[oi].re*twim
//     				     	  + arranged[oi].im*twre;

     				temp.re =  arranged[oi].re*twiddles[tindex].re
     			     		   - arranged[oi].im*twiddles[tindex].im;
     			    temp.im =  arranged[oi].re*twiddles[tindex].im
     			     		   + arranged[oi].im*twiddles[tindex].re;

    				tindex += lasttinc;

     			    arranged[oi].re = arranged[ei].re - temp.re;
     			    arranged[oi++].im = arranged[ei].im - temp.im;
     			    arranged[ei].re += temp.re;
     			    arranged[ei++].im += temp.im;

     			}	// while

     		}	// else

     		// if original data was in-phase only, with no quadrature component,
     		// make a copy of the FFT result to be ready for post-processing
     		// N.B. the underlying array needs to be twice the current array size
     		if(realonly)
     		{
    			int i,j,k,m;

     			// copying originals into second half of underlying array
     			j = start + length;
     			m = start + oddoffset;
     			k = start + oddoffset + length;

     			for(i = start; i < nextstart; i++){
     				arranged[j].re 	 = arranged[i].re;
     				arranged[j++].im = arranged[i].im;

     				arranged[k].re   = arranged[m].re;
     				arranged[k++].im   = arranged[m++].im;
     			}	// for

     		}	// if

     	}	 // else

     	tinc >>= 1;		// as the butterfly size is doubled the steps are halved

    }	// for

}	// fftProc()


// see fft.h for description
void fftMulti(fft_cpx* cin, const int lenP2,
			  const fft_cpx* twiddles,
    		  const int whichcpus, const int normalise,
    		  const int arrange, const int realonly)
{
	/*
	 * Input validation
	 */
	if((FFT_MIN_LENP2 > lenP2)||(FFT_MAX_LENP2 < lenP2))
	{
		printf("\nfftMulti():Invalid power of 2 for array length: %d, should be between %d and %d\n",
				lenP2, FFT_MIN_LENP2, FFT_MAX_LENP2);
		return;
	}	// if

	if ((0 == whichcpus)||(whichcpus & ~CPUS_MASK))
	{
	    printf("\nInvalid choice of CPUs, input : 0x%x, must be between 0 and 0x%x\n", whichcpus, CPUS_MASK);
	    return;
	}	// if


	/*
	 * Processing
	 */
	int length = 1 << lenP2;

	int blocksize;

	int numcpus = getNumCpus(whichcpus);
	int first,second,next;

	fft_cpx * use;
	if(arrange)
	{
		use = (void *)FFT_RAIN;
	} else {
		use = cin;
	}	// else

	switch (numcpus)
	{
		case 1:
			blocksize = 1 << lenP2;			// length

			// main processing
			giveToCpu(whichcpus,
					  cin, lenP2,
					  0, blocksize,
					  0, 1,
					  twiddles, 0, 1 << (lenP2-1),
					  normalise, arrange,
					  realonly,0);

			// MUST WAIT FOR THiS TO HAVE FINISHED
			fftWait(whichcpus);

			if(realonly)
			{
				giveToCpu(whichcpus,
						  use, lenP2,
						  0, blocksize,
						  0, 1,
						  twiddles, 0, 2,
						  normalise, 0,
						  realonly,1);
			}	// if

			// MUST WAIT FOR THiS TO HAVE FINISHED
			fftWait(whichcpus);
			break;

		case 2:
		case 3:
			blocksize = 1 << (lenP2 - 1);	// half of the length

			// get the CPU to use
			if (0 != (whichcpus & 8))first = 8;			// 1xxx
			else if (0 != (whichcpus & 4))first = 4;	// 01xx
			else first = 2;								// 001x

			next = first >> 1;

			if (0 != (whichcpus & next))second = next;
			else if (0 != (whichcpus & (next >> 1)))second = (next >> 1);
			else second = 1;

			// main processing - most of the time is spent here
			giveToCpu(first,
					  cin, lenP2,
					  0, blocksize,
					  0, 1,
					  twiddles, 0, 1 << (lenP2-1),
					  normalise, arrange,
					  realonly, 0);

			giveToCpu(second,
					  cin, lenP2,
					  blocksize, length,
					  0, 1,
					  twiddles, 0, 1 << (lenP2-1),
					  normalise, arrange,
					  realonly,0);

			// MUST WAIT FOR THESE TO HAVE FINISHED
			fftWait(first|second);			// wait for both to finish

			// combine results
			giveToCpu(first,
					  use, lenP2,
					  0, blocksize,
					  blocksize >> 1, lenP2-1,
					  twiddles, 0, 1,
					  normalise, 0,
					  realonly, 0);

			giveToCpu(second,
					  use, lenP2,
					  blocksize >> 1, (blocksize >> 1) + blocksize,
					  blocksize >> 1, lenP2-1,
					  twiddles, blocksize >> 1, 1,
					  normalise, 0,
					  realonly, 0);

			// MUST WAIT FOR THESE TO HAVE FINISHED
			fftWait(first|second);			// wait for both to finish

			// complete the real only processing
			if(realonly)
			{
				giveToCpu(first,
						  use, lenP2,
						  0, blocksize,
						  0, 1,
						  twiddles, 0, 2,
						  normalise, 0,
						  realonly, 1);

				giveToCpu(second,
						  use, lenP2,
						  blocksize, length,
						  0, 1,
						  twiddles, blocksize, 2,
						  normalise, 0,
						  realonly, 1);

				// MUST WAIT FOR THESE TO HAVE FINISHED
				fftWait(first|second);			// wait for both to finish
			}	// if

			break;

		case 4:

			blocksize = 1 << (lenP2 - 2);	// quarter of the length

	    	// main processing - most of the time is spent here
	    	giveToCpu(8,
	    			  cin, lenP2,
	    			  0, blocksize,
	    			  0, 1,
	    			  twiddles, 0, 1 << (lenP2-1),
	    			  normalise, arrange,
	    			  realonly, 0);
	    	giveToCpu(4,
	    			  cin, lenP2,
	    			  blocksize, blocksize << 1,
	    			  0, 1,
	    			  twiddles, 0, 1 << (lenP2-1),
	    			  normalise, arrange,
	    			  realonly, 0);
	    	giveToCpu(2,
	    			  cin, lenP2,
	    			  blocksize<<1, (blocksize << 1) + blocksize,
	    			  0, 1,
	    			  twiddles, 0, 1 << (lenP2-1),
	    			  normalise, arrange,
	    			  realonly, 0);
	    	giveToCpu(1,
	    			  cin, lenP2,
	    			  (blocksize << 1) + blocksize, length,
	    			  0, 1,
	    			  twiddles, 0, 1 << (lenP2-1),
	    			  normalise, arrange,
	    			  realonly, 0);

			// MUST WAIT FOR THESE TO HAVE FINISHED
			fftWait(0xf);			// wait for all to finish

			// combine results
	    	giveToCpu(8,
	    			  use, lenP2,
	    			  0, blocksize,
	    			  blocksize>>1, lenP2 - 2,
	    			  twiddles, 0, 2,
	    			  normalise, 0,
	    			  realonly, 0);
	    	giveToCpu(4,
	    			  use, lenP2,
	    			  blocksize>>1, (blocksize >> 1) + blocksize,
	    			  blocksize >> 1, lenP2 -2,
	    			  twiddles, blocksize, 2,
	    			  normalise, 0,
	    			  realonly, 0);
	    	giveToCpu(2,
	    			  use, lenP2,
	    			  blocksize<<1, 3*blocksize,
	    			  blocksize >> 1, lenP2 -2,
	    			  twiddles, 0, 2,
	    			  normalise, 0,
	    			  realonly, 0);
	    	giveToCpu(1,
	    			  use, lenP2,
	    			  (blocksize <<1) + (blocksize >> 1), 3 * blocksize + (blocksize >> 1),
	    			  blocksize >> 1, lenP2 -2,
	    			  twiddles, blocksize, 2,
	    			  normalise, 0,
	    			  realonly, 0);

			// MUST WAIT FOR THESE TO HAVE FINISHED
			fftWait(0xf);			// wait for all to finish

			// combine again
	        giveToCpu(8,
	        		  use, lenP2,
	        		  0, blocksize,
	        		  3*(blocksize >> 1), lenP2-2,
	        		  twiddles, 0, 1,
	        		  normalise, 0,
	        		  realonly, 0);
	        giveToCpu(4,
	        		  use, lenP2,
	        		  blocksize >> 1,(blocksize >> 1) + blocksize,
	        		  3*(blocksize >> 1), lenP2-2,
	        		  twiddles, blocksize>>1, 1,
	        		  normalise, 0,
	        		  realonly, 0);
	        giveToCpu(2,
	        		  use, lenP2,
	        		  blocksize,(blocksize<<1),
	        		  3*(blocksize >> 1), lenP2-2,
	        		  twiddles, blocksize, 1,
	        		  normalise, 0,
	        		  realonly, 0);
	        giveToCpu(1,
	        		  use, lenP2 ,
	        		  blocksize + (blocksize >> 1),(blocksize << 1) + (blocksize >>1),
	        		  3*(blocksize >> 1), lenP2-2,
	        		  twiddles, 3*(blocksize >> 1), 1,
	        		  normalise, 0,
	        		  realonly, 0);

			// MUST WAIT FOR THESE TO HAVE FINISHED
			fftWait(0xf);			// wait for all to finish

			// complete the real only processing
			if(realonly)
			{
				giveToCpu(8,
						  use, lenP2,
					      0, blocksize,
					      0, 1,
					      twiddles, 0, 2,
					      normalise, 0,
					      realonly, 1);

				giveToCpu(4,
						  use, lenP2,
					      blocksize, blocksize << 1,
					      0, 1,
					      twiddles, 0, 2,
					      normalise, 0,
					      realonly, 1);

				giveToCpu(2,
						  use, lenP2,
					      blocksize<<1, 3*blocksize,
					      0, 1,
					      twiddles, 0, 2,
					      normalise, 0,
					      realonly, 1);

				giveToCpu(1,
						  use, lenP2 ,
					      3*blocksize,blocksize << 2,
					      0, 1,
					      twiddles, 0, 2,
					      normalise, 0,
					      realonly, 1);

				// MUST WAIT FOR THESE TO HAVE FINISHED
				fftWait(first|second);			// wait for both to finish
			}	// if

			break;

	}	// switch

}	// fftMulti()


// see fft.h for description
void giveToCpu(const int which,								// 1,2,4 or 8
			   fft_cpx* cin, const int lenP2,				// input array
			   const int begin, const int nextbegin,		// part indices
			   const int oddadd, const int startP2,
			   const fft_cpx* twiddles,
			   const int toffset, const int twinc,
			   const int normalise, const int arrange,
			   const int realonly, const int secondpass)
{
	if(0 == which) return;					// nothing to do this time


	if((FFT_MIN_LENP2 > lenP2)||(FFT_MAX_LENP2 < lenP2))
	{
		printf("\ngiveToCpu(): Invalid power of 2 for array length: %d, should be from %d to %d\n",
				lenP2, FFT_MIN_LENP2, FFT_MAX_LENP2);
		return;
	}	// if

	// get the index into the FFT_FLAGS
    int cpuindex;
    if(1 == which) cpuindex = 0;
    else if(2 == which) cpuindex = 1;
    else if(4 == which) cpuindex = 2;
    else if(8 == which) cpuindex = 3;

    else return;

	volatile fft_cpu *thisCpu = (void *)(FFT_FLAGS + cpuindex*sizeof(fft_cpu));

	thisCpu->cinadrs	= (unsigned int)cin;
	thisCpu->lenP2		= (unsigned int)lenP2;
	thisCpu->begin		= (unsigned int)begin;
	thisCpu->nextbegin	= (unsigned int)nextbegin;
	thisCpu->oddadd		= (unsigned int)oddadd;
	thisCpu->startP2	= (unsigned int)startP2;
	thisCpu->twadrs 	= (unsigned int)twiddles;
	thisCpu->toffset	= (unsigned int)toffset;
	thisCpu->tinc		= (unsigned int)twinc;
	thisCpu->normalise 	= (unsigned int)normalise;
	thisCpu->arrange 	= (unsigned int)arrange;
	thisCpu->realonly   = (unsigned int)realonly;
	thisCpu->secondpass = (unsigned int)secondpass;

	thisCpu->spare2		= 0xa5a5a5a5;					// for debug
	thisCpu->status		= (unsigned int)FFT_START;		// do this last

    // if the CPU is the one running this code, do the FFT here
    if(START_CPU_INDEX == cpuindex)
    {
    	thisCpu->status |= FFT_START;
		thisCpu->spare1 = 0xbabecafe;				// for debug

    	fftProc((void *)thisCpu->cinadrs, thisCpu->lenP2,
				thisCpu->begin, thisCpu->nextbegin,
				thisCpu->oddadd, thisCpu->startP2,
				(void *)thisCpu->twadrs,
				thisCpu->toffset, thisCpu->tinc,
				thisCpu->normalise, thisCpu->arrange,
				thisCpu->realonly, thisCpu->secondpass);

    	thisCpu->spare2 = 0xcafebabe;				// for debug
		thisCpu->status &= ~FFT_WORK_MASK;			// finished

    } else {
    	*((int *)IRQMP_ADRS) = 0x380b0000 | which;
    }	// else

}	// giveToCpu()


// see fft.h for description
void fftWait(int whichcpus){

	// other CPUS loop forever, check for work, do it, set done status

	volatile fft_cpu *cpu0 = (void *)(FFT_FLAGS);
	volatile fft_cpu *cpu1 = (void *)(FFT_FLAGS + sizeof(fft_cpu));
	volatile fft_cpu *cpu2 = (void *)(FFT_FLAGS + 2*sizeof(fft_cpu));
	volatile fft_cpu *cpu3 = (void *)(FFT_FLAGS + 3*sizeof(fft_cpu));

//	volatile fft_cpu * cpus = (void *)FFT_FLAGS;
	int notdone = 0;
	do
	{
		notdone =  (8 & whichcpus)? FFT_WORK_MASK & cpu3->status : 0;
		notdone |= (4 & whichcpus)? FFT_WORK_MASK & cpu2->status : 0;
		notdone |= (2 & whichcpus)? FFT_WORK_MASK & cpu1->status: 0;
		notdone |= (1 & whichcpus)? FFT_WORK_MASK & cpu0->status: 0;

	} while (notdone);

	return;

}	// fftWait()

// see fft.h for description
void fftRealReturn(fft_cpx * cin, const int halflenP2, fft_cpx * twiddles)
{
	fft_scalar value;
	int i,j;

	int length = 1 << halflenP2;

	// results overwrite originals that are needed subsequently,
	// so copy originals into second half of underlying array
	j = length;
	for(i =0; i <= length; i++){
		cin[j].re 	 = cin[i].re;
		cin[j++].im  = cin[i].im;
	}	// for

	j = (length << 1) +1;
	for(i = 0; i < length; i++)
	{
		j--;
		value =  (cin[i].re + cin[j].re
			      + twiddles[i].re*(cin[i].im + cin[j].im)
			      + twiddles[i].im*(cin[i].re - cin[j].re))/2.0;

		cin[i].im = (cin[i].im - cin[j].im
			         + twiddles[i].im*(cin[i].im + cin[j].im)
			         - twiddles[i].re*(cin[i].re - cin[j].re))/2.0;

		cin[i].re = value;

	}	// for

	// Nyquist case
	// if they existed, twiddles[length].re would be -1, twiddles[length].im would be 0;
	cin[length].re = cin[length].re - cin[length].im;
	cin[length].im = 0;

	// copy positive frequency  results to negative, doing complex conjugate
	for(i = 1; i < length ; i++)
 	{
 		cin[length+i].re = cin[length-i].re;
 		cin[length+i].im = -cin[length-i].im;
 	}	// for

	return;
}	// fftRealReturn()


/*
 * UTILITIES
 */

// see fft.h for description
int getNumCpus(const int whichcpus)
{
	// set up count of cpus to use
	int numcpus = 0;				// how many cpus to use

	numcpus += (whichcpus & 1)? 1 : 0;
	numcpus += (whichcpus & 2)? 1 : 0;
	numcpus += (whichcpus & 4)? 1 : 0;
	numcpus += (whichcpus & 8)? 1 : 0;

    return numcpus;

}	// getNumCpus


// see fft.h for description
int settingsIncorrect(const int lenP2, int *cpusrequest, const int numruns, const int type, int* lines)
{
	int result = 0;
	// array size is limited by the way memory is being allocated at present
	if((FFT_MIN_LENP2 > lenP2)||(FFT_MAX_LENP2 < lenP2))
	{
		printf("\nInvalid power of 2 for array length: %d, should be between %d and %d\n",
				lenP2, FFT_MIN_LENP2, FFT_MAX_LENP2);
		result |= 1;
	}	// if

	int whichcpus = *cpusrequest;
	if((0 == whichcpus)||(whichcpus & ~CPUS_MASK))
	{
		printf("\nInvalid choice of CPUs: 0x%x, must be from 1 to 0x%x\n",
				whichcpus, CPUS_MASK);
		result |= 1;
	}	// if

	// get the count of cpus being requested
	int numcpus = getNumCpus(whichcpus);

	// only allow 1, 2 or 4 cpus, if 3 deselect the lowest order one
	if(3 == numcpus)
	{
	    if (whichcpus & 1) whichcpus &= 0xe;		// deselect CPU 0
	    else whichcpus &= 0xd;						// deselect CPU 1
	    *cpusrequest = whichcpus;
	    numcpus = 2;
	    printf("\nChoice of 3 CPUs invalid, will use these 2 CPUs : 0x%x\n", whichcpus);
	}	// if

	// find how many CPUs are present on system
	volatile int status = *(int *)IRQMP_ADRS;		// current system configuration
	int statuscpus = (((status & CPUS_PM1_MASK)>>CPUS_PM1_SHIFTS)&0xf) + 1;
	printf("\nMain, Status: 0x%x, CPUs present %d,  CPUs requested: %d, WhichCPUs: 0x%x\n",
			status,statuscpus,numcpus,whichcpus);	// debug

	if((0 == numcpus)||(numcpus > statuscpus))		// sanity check
	{
		printf("Abandoning, number of CPUs requested: %d, number present %d",numcpus, statuscpus);
		result |= 1;
	}	// if

	// check the type
	if((0 != type)&&(1 != type))
	{
		printf("Abandoning, invalid type input: %d, should be 0 for forward FFT or 1 for inverse FFT\n",type);
		result |= 1;
	}	// if

	// check the choice of number of lines to display
	if(0 >= *lines)
	{
		printf("Number of lines chosen was %d, no data or result will be displayed\n",*lines);
	} else if ((1 << lenP2) < *lines){
		printf("Number of lines chosen %d exceeds array length %d, array length will be used\n",*lines,1<<lenP2);
		*lines = 1 << lenP2;
	}
	return result;

}	// settingsIncorrect()


// see fft.h for description
void showC(const fft_cpx* in, const int length, const double tolerance)
{
	int i=0;
	float temp1 = 0.0;
	float temp2 = 0.0;
	printf("Complex numbers array, showing first %d elements\n", length);

	for (i=0; i<length; i++) {

	    temp1 = in[i].re;
	    if(  ((0.0 < temp1)&&( tolerance > temp1))
	       ||((0.0 > temp1)&&(-tolerance < temp1))){
	   	   temp1 = 0.0;
	    }	// if

	    temp2 = in[i].im;
	    if(  ((0 < temp2)&&( tolerance > temp2))
	       ||((0 > temp2)&&(-tolerance < temp2))){
	  	   temp2 = 0.0;
	    }	// if

	    printf("%d:   %g      %g\n", i, temp1, temp2);
	}	// for

}	// showC()


// see fft.h for description
void show8C(const fft_cpx* in, const int length, const double tolerance, const int realonly)
{
	if(32 > length)
	{
		showC(in,length,tolerance);
	}	// if

	int i,k;
	fft_cpx temp;
	float temp1 = 0.0;
	float temp2 = 0.0;
	printf("Complex numbers array, showing certain elements\n\n");

	int nyindex = length>>1;
	int begins[4] = {0, nyindex - 8, nyindex, length-8};

	for(k =0; k<4; k++)
	{
		for (i=0; i<8; i++) {

			temp  = in[begins[k]+i];
			temp1 = temp.re;
			if(  ((0.0 < temp1)&&( tolerance > temp1))
			   ||((0.0 > temp1)&&(-tolerance < temp1))){
		   	   temp1 = 0.0;
			}	// if

			temp2 = temp.im;
			if(  ((0 < temp2)&&( tolerance > temp2))
		       ||((0 > temp2)&&(-tolerance < temp2))){
		  	   temp2 = 0.0;
		    }	// if

		    printf("%d:\t\t\t%g      %g\n", begins[k]+i, temp1, temp2);

		    if((0 == begins[k]+i)||(nyindex == begins[k]+i)) printf("\n");

		    if(realonly && (nyindex == begins[k] + i))
		    {
		    	printf("In-phase data only, last frequency shown above is Nyquist\n");
		    	printf("Negative frequencies not shown (complex conjugates of positive frequencies)\n\n");
		    	break;	// for loop
		    }	// if

		}	// for

		printf("\n");
		if(realonly && (2 == k)) break;

	}	// for

}	// show8C()
