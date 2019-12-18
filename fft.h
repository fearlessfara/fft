/*
 ============================================================================
 Name        : fft.h
 Author      : Faraone Christian Gennaro & Professor Michael Ryan
 Version     : 1.0
 Copyright   : Copyright OCE Technology
 Description : Multicore FFT

 ============================================================================
 */
#ifndef FFT_H_
#define FFT_H_

#include <math.h>
#include <time.h>
#include <unistd.h>

// assuming basic number type is float (for now)
# ifndef fft_scalar
#   define fft_scalar float
# endif

// (in-phase, quadrature) data is treated as a complex number
typedef struct {
    fft_scalar re;					// real, i.e. in-phase
    fft_scalar im;					// imaginary, i.e. quadrature
}fft_cpx;

// each CPU has a struct of this type set up in SRAM
typedef struct {
    unsigned int status;
    unsigned int cinadrs;
    unsigned int lenP2;				// power of 2 that gives length
    unsigned int begin;
    unsigned int nextbegin;
    unsigned int oddadd;
    unsigned int startP2;			// power of 2 that gives butterfly to start with
    unsigned int twadrs;
    unsigned int toffset;
    unsigned int tinc;
    unsigned int normalise;
    unsigned int arrange;
    unsigned int realonly;
    unsigned int secondpass;
    unsigned int spare1;			// bring up to a 16 word boundary
    unsigned int spare2;
}fft_cpu;

#define FFT_WORK_MASK   1			// pick off last bit
#define FFT_NOWORK		0			// nothing to be done/finished
#define FFT_START		1			// work to be done/working

#define MAX_CPUS 		4
#define LAST_CPU		8			// 0x1000
#define START_CPU_INDEX	0			// N.B. fft_multi() assumes is 0, CPU that starts on power up

#define IRQMP_ADRS		0x80000210  // address of multiprocessor status register
#define CPUS_PM1_MASK	0xf0000000	// tells how many CPUs in system minus 1
#define CPUS_PM1_SHIFTS 28			// shift them to low order
#define	CPUS_MASK    	0xf			// bits used to power-up CPUs
#define CPU_START_PAUSE	1000		// pause between and after CPU startups

// set up space for the various arrays. S698PM SRAM is 32 Mbyte from 0x40000000.
#define FFT_SRAM 	 	0x60000000				// half way
#define FFT_FLAGS	 	FFT_SRAM+0x80000
#define FFT_CIN  	 	FFT_SRAM+0x100000		// space for 2 power 20 bytes,
#define FFT_RAIN	 	FFT_SRAM+0x200000		// i.e. 2 power 18 doubles
#define FFT_TWIDDLES 	FFT_SRAM+0x400000		// i.e. 2 power 17 complex nos
#define FFT_COUT	 	FFT_SRAM+0x600000
#define FFT_MAP			FFT_SRAM+0x800000		// maps input index to rearranged index

#define FFT_MIN_LENP2	5			// minimum input 32 (in-phase,quadrature)
#define FFT_MAX_LENP2	17			// for consistency with space allocation above

/*
 * FUNCTION DECLARATIONS
 */

/*
 *  The functions declared below are used to create an FFT
 *  using 1, 2 or 4 processors acting in parallel.
 *
 *	At present the input array length must be a power of 2.
 *
 *  At present the input data in each (in-phase, quadrature) pair must be C floats.
 *
 *	The processing is based on arrays of complex numbers corresponding to
 *	the (in-phase, quadrature) data.
 *
 *	If the input data is in-phase only, the input array will have 0 in each entry's
 *	quadrature component. An option is provided to allow this be exploited so as
 *	to get the result roughly twice as fast
 *
 *  Usage:
 *		0.	Create an array holding the in-phase and quadrature data
 *			as complex numbers (fft_cpx).
 *
 *			At present, the length of this array must be a power of 2.
 *			At present, data elements should be of the C float type.
 *
 *		1.	The complex array must be arranged in an order suitable
 *			for 'discrimination in time' (DIT) processing.
 *
 *			This can be done as the entries are made one at a time.
 *			An index array can first be created for use to map
 *			successive input to the correct position in the
 *			array. This can be stored in ROM, and reused to allocate data
 *			as it comes in to the appropriate place in the array.
 *
 *			The function fftArrangeIndex() below is provided for
 *			this purpose.
 *
 *			If the input data is not rearranged as it arrives,
 *			a parameter can be set to cause rearrangement to be done
 *			in parallel in the first phase of FFT processing on each CPU.
 *			In that case the FFT result does not overwrite the input data,
 *			but can be found in the FFT_RAIN area defined below.
 *
 *			Alternatively, if the array already exists (as an array of complex numbers),
 *			it can be rearranged in situ using fftArrangeIn(),
 *			or as a new complex array using fftArrangOut().
 *
 *			The data in the rearranged array is overwritten by the FFT
 *			results, so if it is desired to keep the original data
 *			either the parallel rearrangement parameter should be set
 *			or the fftArrangeOut() function should be used.
 *
 *			To avoid the overhead of rearrangement it is best to rearrange
 *			the data as it arrives, if possible.
 *
 *		3.	The 'twiddles' array, containing complex roots of unity,
 *			must be created. This is half the length of the input
 *			array of complex numbers, as due to symmetry only half
 *			of the roots need to be created and stored.
 *			The function fillTwiddles() is provided for this purpose.
 *
 *		N.B.	If the size of the data array does not change,
 *				the twiddles array can be stored in ROM.
 *				Two versions are needed in that case, one for the forward FFT,
 *				one for the inverse FFT. Both can be generated using fillTwiddles().
 *
 *		4.	The appropriately arranged complex input array is passed to fftMulti(),
 *			together with the twiddles array and which processors to use.
 *			The result overwrites the original data in the rearranged array.
 *
 *		5.	Other than input rearrangement, no copying is involved in passing data
 *			to the different processors, each is given a reference to a contiguous
 *			portion of the input array.
 *			A reference to the twiddles array is also passed, though the entries
 *			used in this by a processor are not contiguous.
 *
 *		6.	All stages of processing up to the last two (4 CPUs) or last one (2 CPUs)
 *			are done independently on the CPUs, with no interaction between them.
 *			The last stage or stages require results generated by a different CPU,
 *			so must wait for the previous stages to complete before progressing.
 *
 *		7.	If the 'real only', i.e in-phase only, option is selected any 'imaginary',
 *			i.e quadrature, component in the input is ignored and the input array compressed
 *			into half its length with the entries in the remaining half set to (0.0,0.0).
 *			The processing is done roughly twice as fast.
 *			Only the positive frequency components are produced,
 *			the negative frequency components are easily derived if needed by taking
 *			the complex conjugates of the positive frequency components,
 *			i.e. changing the sign of the imaginary component and leave the rest the same.
 */

/*
 * FUNCTION DECLARATIONS
 */
/*
 * This ASM-function returns the CPU index (0 to 3) of the current S698PM CPU.
 */
 static inline unsigned int get_asr17(void);


/*
 * Starting CPUs
 *
 * All CPUs have the same entry point.
 * The linker defaults to 0x40000000.
 * DMON's ep command is used to set each CPU's initial PC to this.
 *
 * This function is entered immediately on entering main().
 * It determines the index (0 to 3) of the CPU running the code.
 * If this is the START_CPU_INDEX it initializes the data areas used by all CPUs
 * and returns to main().
 * Otherwise, it checks to see if any work is waiting, does it,
 * and halts the CPU. When work becomes available, the CPU is restarted,
 * does the work, and halts again. This function never returns if running on
 * these CPUs.
 *
 */
int fftCpusStart(void);


/*
 * Preliminary start of the selected CPUs if other than the one running this.
 * Each such CPU starts, runs through its initialization code, and powers down
 * Each CPU shares the same initialization code, so there is a pause between
 * each startup just in case ...
 *
 * A CPU is restarted when there is work for it to do
 * (restarts are done directly, not by this).
 *
 *  Inputs:
 *  	const int		which cpus		from 0 to 0xf, the current CPU, if selected, is ignored
 *
 *  Returns:							None
 *
 * 						N.B. This routine should not be used to restart CPUs
 */
void startCpus(const int whichcpus);


/*
 * INPUT ARRAY REARRANGEMENT
 *
 * This is used with a discrimination in time (DIT) implementation of the FFT.
 * Separate parts of the rearranged array can be processed independently
 * on different CPUs without any interaction between the CPUs being required.
 * The last two stages (4 CPUS) or just the last stage (2 CPUS)can only be done when
 * the previous stages are finished, and are again distributed across all CPUs.
 *
 * The best approach to input rearrangement is to do it as the data arrives.
 *
 * Alternatively, by setting a parameter rearrangement can be done in parallel
 * one each CPU.
 *
 * Three very similar functions, not parallel, are given to allow a complete input
 * array be rearranged:
 * 	fftArrangeOut		an out of place rearrange, the input array is unchanged
 * 	fftArrangeIn		an in place rearrange, the input array is rearranged in situ
 * 	fftaArrangeIndex	creates an array of indexes. Entry [i] gives the index where
 * 						the (in-phase, quadrature) complex pair should be stored when
 * 						data is read one pair at a time
 *						- the input array will then not need further rearrangement
 *
 * All three functions assume the array length is a power of 2, and this power is what
 * is input to give the array length.
 *
 * (Thanks for bit-swap method used to graphics.stanford.edu/~seander/bithacks.html#ReverseParallel)
 *
 */

/*
 *	Out of place rearrangement of input array to prepare for FFT.
 *	Assumes array length is a power of 2.
 *
 *	Inputs:
 *		fft_cpx* 		cin		pointer to an array of complex numbers
 *								(the length of this array must be a power of 2)
 *		int				lenP2	log to base 2 of array length. 0 <= lenP2 <= 30.
 *
 *		fft_cpx* 		result	pointer to resulting array of complex numbers
 *								(same length as cin array)
 *
 *	Result:				A rearranged version of cin input array in the result array
 *
 *	N.B.				Space for the result array of complex numbers must be allocated
 *						before calling this function, same length as the cin input array.
 *
 */
void fftArrangeOut(fft_cpx* cin, fft_cpx* result, int lenP2);


 /*
  *	In place rearrangement of input array to prepare for FFT.
  *	Assumes array length is a power of 2.
  *
  *	Inputs:
  *		fft_cpx* 		cin		pointer to an array of complex numbers
  *								(the length of this array must be a power of 2)
  *		int				lenP2	log to base 2 of array length.
  *
  *	Result:				A rearranged version of cin input array in the same space
  *
  */
 void fftArrangeIn(fft_cpx* x, int lenP2);


 /*
  *	Generate an array of indices giving the array position where each successive input
  *	should be placed in the input array.
  *	This array could be stored in ROM and used when inputting data so as to avoid
  *	any subsequent need to spend time rearranging the input array.
  *
  *	Assumes amount of input datais a power of 2.
  *
  *	Inputs:
  *		int * 			mapindex	pointer to an array of ints that will hold the index to be used
  *									in the complex input array for each successive data input.
  *									(the length of this array must be a power of 2
  *		int				lenP2		log to base 2 of array length. 0 <= lenP2 <= 30
  *									(regarding the array as an array of complex numbers)
  *
  *	Result:				mapindex will hold appropriate index in data array to use to hold complex input i
  *
  *	N.B. if the input data is an array of scalars, the length of mapindex[] will be half of the length
  *		 of the array to hold the scalars, and the input index k used to look up the mapindex array should only
  *	     increment after each pair of reals have been added to dataarray[mapindex[k]<< 1] and dataarray[(mapindex[k]<< 1) +1]
  *	     (they will be treated as the real and imaginary parts of the complex number at the index given by mapindex[k])
  *
  *	N.B. Space for the mapindex array of ints must be allocated before calling this function,
  *		 same length as complex data input array, half length of real data input arrays.
  *
  *	N.B. The mapping array is its own inverse
  *
  */
 void fftArrangeIndex(int* mapindex, int lenP2);


 /*
   *	If the 'realonly' option is selected, the input data is treated as having in-phase components only
   *	and any quadrature component is ignored and will be overwritten in processing.
   *	The 'real', i.e. in-phase input components are rearranged into half the number of inputs
   *	with every second in-phase input now in the quadrature position, and zero in remaining half of the input array.
   *	This allows faster FFT processing but requires post-processing.
   *	(Post processing is distributed across all processors and done later in a final processing phase.
   *	 It produces the positive frequency components only, the negative components can be obtained by taking
   *	 the complex conjugates of these.)
   *
   *	Assumes amount of input data is a power of 2.
   *
   *	Inputs:
   *		fft_cpx * 		cin			input array with all 0s in quadrature (imaginary) positions.
   *									Altered so that quadrature positions hold alternate real inputs,
   *									and remaining array entries are made (0.0,0.0)
   *		int				lenP2		log to base 2 of array length. 0 <= lenP2 <= 30.
   *									(regarding the array as an array of complex numbers)
   *
   *	Result:				first half of cin holds rearranged data, rest is (0.0,0.0)
   *
   *						N.B. 		use of this function can usually be avoided by simply casting an input
   *									array of scalars to be an array of complex numbers.
   */
 void fftRealOnlyArrange(fft_cpx * cin, const int lenP2);


 /*
  * TWIDDLES ARRAY
  */

 /*
  * For an input data array of length N the 'twiddle' factors are the N complex roots
  * of unity. Only N/2 of these are calculated, as the others are just minus these,
  * and are provided by using subtraction rather than addition.
  *
  * Getting the complex roots of unity is time consuming, so where the size of
  * array to be processed is known in advance the twiddles array should be stored
  * in ROM. Two arrays should be stored, one for the forward FFT, one for the inverse.
  * Each has length N/2, half that of the complex data array.
  *
  * Time consuming calculations of sines and cosines are involved, but in some cases
  * can be avoided by exploiting various symmetries. These require the complex data array
  * size to be at least 8, but it is restricted here to at least 1 << FFT_MIN_LENP2, or 32 at present.
  *
  *
  * Inputs:
  *		fft_cpx* 		twiddles	pointer to an array of complex numbers
  *									(the length of this array will be half the length
  *									 of the complex data array with which the twiddles
  *									 array will be used)
  *		const int		lenP2		log to base 2 of complex data array length.
  *									(from FFT_MIN_LENP2	to FFT_MAX_LENP2, i.e. 5 to	17)
  *
  *		const int		type 		0 for FFT, 1 for inverse FFT
  *
  * Result:				N/2 of the N complex roots of unity in twiddles array
  *
  * N.B.				Space for the twiddles array of N/2 complex numbers must be allocated
  *						before calling this function, where N is the length of the complex data
  *						array. i.e. 1 << lenP2.
  *
  */
 void fillTwiddles(fft_cpx* twiddles, const int lenP2, const int type);


 /*
  * FFT PROCESSING
  *
  * The (in-phase, quadrature) input data is regarded as an array of complex numbers.
  * The resulting FFT is also an array of complex numbers.
  *
  * If the data is not already in rearranged order, by setting the 'arrange' parameter
  * the input array can be rearranged in parallel, with each CPU doing its own section.
  * This involves copying the original input data into the area at FFT_RAIN,
  * where it will eventually be overwritten by the FFT result, with the original
  * input data left unchanged.
  *
  */

 /*
  * This function calculates all or part of the FFT.
  * It can be used on separate processors or in separate threads.
  *
  * It processes part or all of a complex data array, replacing parts of the
  * array with appropriate FFT data. After any initial rearrangement
  * all processing is done in situ .
  *
  * The length of a sub-array must be a power of 2 and no less than 8.
  * The sub-array can be the whole array.
  *
  * The size of the butterfly blocks with which to start processing a sub-array is
  * input as a power of 2.
  *
  * Where the blocks in the subarray have not already been processed,
  * the starting power of 2 is 1.
  *
  * If parts of a subarray have already been processed, the correct butterfly
  * step must be given. All parts of the subarray must have the same size.
  *
  * The post-processing phrase needed when the original data is in-phase only
  * is done if 'secondpass' is set.
  *
  * Rearrangement of the data for this section is done if the 'arrange' is set.
  *
  *
  * Inputs:
  *		fft_cpx* 		arranged	pointer to the array of complex numbers that
  *									initially contains the rearranged complex data input.
  *									Will hold result.
  *		const int		lenP2		Power of 2 that gives length of input array.
  *
  *		const int		start		position in array at which sub-array starts
  *		const int		nextstart	position one after end of data to process
  *									N.B. (nextstart - start) must be a power of 2
  *
  *		const int		oddadd		Added internally to half of sub-array size.
  *									Used to get distance from even to odd entries.
  *
  *		const int		startP2		power of 2 that gives butterfly at which to start.
  *									Must be at least 1.
  *									1 if no previous processing done on this sub-array.
  *									Otherwise gives size of butterfly to start with.
  *
  *		const fft_cpx*	twiddles	pointer to twiddles array - this must correspond
  *									to the original complete complex input array.
  *		const int		toffset		initial offset into twiddles array
  *		const int		twinc		initial increment to using going through twiddles array
  *
  *		const int		normalise	0 for don't normalize, anything else -> normalize
  *
  *		const int		arrange		0 for don't arrange, anything else -> arrange
  *
  *		const int		realonly	0 if original data complex, 1 if in-phase only
  *									and has been rearranged
  *
  *		const int		secondpass	0 if not ready to do post-processing,
  *									1 to do post-processing only
  *
  * Result:				FFT processed values in situ in this subarray
  *
  *	N.B.				If 'realonly' is selected, assumes that the input array is a compressed
  *						version of something twice as long, and that this additional space
  *						is available.
  *
  *	N.B.				If 'realonly' is selected, FFT result contains the positive frequencies
  *						only. Negative frequencies are easily derived if needed as the complex
  *						conjugates of these.
  *
  */
 void fftProc(fft_cpx* arranged, const int lenP2,
 			 const int start, const int nextstart,
 			 const int oddadd, const int startP2,
     		 const fft_cpx* twiddles,
     		 const int toffset, const int twinc,
     		 const int normalise, const int arrange,
     		 const int realonly, const int secondpass);


 /*
  * Allocates work of doing FFT to one or more CPUs
  * and combines results from them to give FFT.
  *
  * The input is an array of complex numbers representing the
  * (in-phase,quadrature) data.
  * Its length must be a power of 2 - this power rather than the
  * length is input.
  *
  * The complex number array giving the FFT is returned in situ
  * in the input array, unless the 'arrange parameter is set,
  * in which case the result is in the area at FFT_RAIN.
  *
  * At present the number of CPUs to use must be 1, 2 or 4.
  *
  * Inputs:
  *		fft_cpx* 		arranged	pointer to an array of complex numbers that
  *									holds the rearranged complex data input.
  *									Will hold result.
  *		const int		lenP2		power of 2 that gives size of array
  *
  *		const fft_cpx*	twiddles	pointer to twiddles array
  *									- this must correspond to input array
  *									  (it will be half its size)
  *
  *		const int		whichcpus	1 to 0xf, which CPUs to use (binary xxxx, x=1 for use)
  *
  *		const int		normalise	0 for don't normalize, anything else -> normalize
  *
  *		const int		arrange		0 for don't arrange, anything else -> arrange
  *
  *		const int		realonly	0 if original data has both components, 1 if in-phase only
  *									and has been rearranged
  *
  * Result:				Complex valued FFT in situ in arranged
  *
  * N.B.				Assumes START_CPU_INDEX is 0
  *
  */
 void fftMulti(fft_cpx* cin, const int lenP2,
 			  const fft_cpx* twiddles,
     		  const int whichcpus,
     		  const int normalise,
     		  const int arrange,
     		  const int realonly);


 /*
  * Gives one part of the FFT calculation to the selected CPU for execution.
  *
  * The input is part of an array of complex numbers.
  * Its length must be a power of 2 - this power rather than the
  * length is input.
  *
  * The output is placed in the appropriate locations in the input array,
  * unless the 'arrange' parameter is set, in which case the output is
  * in the area at FFT_RAIN.
  *
  * At present the identifier of the CPU to use must be 1, 2, 4 or 8, anything
  * else the function returns without doing anything.
  *
  * Inputs:
  *		const int		which		which CPU to use - must be 1, 2, 4 or 8
  *
  *		fft_cpx* 		arranged	pointer to an array of complex numbers that
  *									holds the rearranged complex data input.
  *									Will hold result.
  *		const int		lenP2		power of 2 that gives size of array
  *
  *		const int		begin		index into arranged, location where this processing is to start
  *
  *		const int		nextbegin	index into arranged, begin <= index < nextbegin
  *
  *		const int		oddadd		Added internally to half of sub-array size.
  *									Used to get distance from even to odd entries.
  * 	const int		startP2		1 << (lenP2 - startP2) gives initial step size in twiddles array
  *
  *		const fft_cpx*	twiddles	pointer to twiddles array
  *									- this must correspond to input array
  *									  (it will be half its size)
  *
  *  	const int		toffset		initial offset into twiddles array
  *		const int		twinc		initial increment to using going through twiddles array
  *
  *		const int		normalise	0 for don't normalize, anything else -> normalize
  *
  *		const int		arrange		0 for don't arrange, anything else -> arrange
  *
  * 	const int		realonly	0 if original data complex, 1 if in-phase only
  *									and has been rearranged
  *
  * Result:				Updated complex values for this stage of FFT processing in arranged
  *
  */
 void giveToCpu(const int which,							// 1,2,4 or 8
 			   fft_cpx* cin, const int lenP2,				// input array
 			   const int begin, const int nextbegin,		// part indices
 			   const int oddadd, const int startP2,
 			   const fft_cpx* twiddles,
 			   const int toffset, const int twinc,
 			   const int normalise, const int arrange,
 			   const int realonly, const int secondpass);


 /*
  * Wait for selected CPUs to have finished.
  * Returns when no CPU has work pending.
  * Can wait indefinitely. TODO put in a time check
  *
  * All stages of the FFT except at most the last two (for the S698PM) proceed independently
  * with no need for communication between CPUs.
  *
  * If using 4 CPUs, each of the last two stages must wait for all CPUs to have finished
  * the preceding stage.
  *
  * If using 2 CPUs, they must wait to process the last stage until both CPUs have finished
  * the earlier stages.
  *
  * Inputs:
  *		const int		whichcpus	which CPUs to check - xxxx, 1 to check, 0 not to check
  *
  *
  * Result:				Wait, only return when no selected CPU has work pending
  */
 void fftWait();


/*
 * UTILITIES
 */

 /*
  * Get the number of CPUs enabled by the selections in whichcpus.
  *
  * 	Input:
  * 		const int 	whichcpus	0 to 0xf  (binary xxxx, 1 = use CPU
  *
  * 	Output:
  * 		int			Count of currently selected CPUs ( 0 to 4)
  *
  */
 int getNumCpus(const int whichcpus);


 /*
  * Validate the input settings and correct if possible
  *
  *	Inputs:
  *		const int		lenP2			must be from FFT_MIN_LENP2	to FFT_MAX_LENP2, (i.e. 5 to 17)
  *
  *		int *			cpusrequest		should be from 0 to 0xf, corresponding to 1, 2 or 4 cpus (but not 3)
  *										(correction attempt made if three CPUs chosen)
  *
  *		const int		numruns			no validation
  *
  *		const int		type			must be 0 or 1
  *
  *		int *			lines			corrected to 1 << lenP2 if greater than this
  *
  *	Returns:
  *		int				0 if all o.k., 1 if a problem (messages given),
  *
  * Inputs:
  *
  */
 int settingsIncorrect(const int lenP2, int *cpusrequest, const int numruns, const int type, int* lines);


 /*
  * Display an array of complex numbers on the console.
  * Can be used to display data input or FFT result.
  *
  * N.B. At present assumes floats being used
  *
  * Inputs:
  * 	fft_cpx* 		in			input array of complex numbers
  *
  *		const int		length		length of the array
  *
  *		const double	tolerance	absolute values less than this are shown as 0.0
  *
  *	Result:				Array is output to the console
  *
  */
 void showC(const fft_cpx* in, const int length, const double tolerance);


 /*
   * Display main parts of FFT result on the console.
   * Shows initial 8 positive frequency components including 0 component,
   * show 8 components before Nyquist frequency, 8 components from Nyquist,
   * and last 8 components.
   * If 'realonly' is selected, shows result up to Nyquist frequency only.
   *
   * N.B. At present assumes floats being used
   *
   * Inputs:
   * 		fft_cpx* 		in			input array of complex numbers
   *
   *		const int		length		length of the array
   *
   *		const double	tolerance	absolute values less than this are shown as 0.0
   *
   *		const int 		realonly	0 to show full FFT, 1 to show non-negative frequencies only.
   *
   *	Result:				Array is output to the console
   *
   */
 void show8C(const fft_cpx* in, const int length, const double tolerance, const int realonly);


#endif // FFT_H_
