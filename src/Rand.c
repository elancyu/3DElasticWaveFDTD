/* Random number generator: MT19937 */
// Copy from MC src directly
// Call: RandR() to generate a random number in uniformly distributed in (0,1)
// Call: RandG() to generate a random number in Gaussian distribution.

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

/* The array for the state vector */
static unsigned long long mt[NN];
/* mti==NN+1 means mt[NN] is not initialized */
static int mti=NN+1;

/* initializes mt[NN] with a seed */
void InitRand(unsigned long long seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++)
        mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}


// generate a random number(real) between 0 and 1
// in: void
// return: the generated number
double RandR(void)
{
    int i;
    unsigned long long x;
    static unsigned long long mag01[2]={0ULL, MATRIX_A};

    if (mti >= NN) { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN+1)
            InitRand(5489ULL);

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

        mti = 0;
    }

    x = mt[mti++];

    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);

    return (x>> 11) * (1.0/9007199254740991.0);
}

// Generate random number following gaussian distribution using Box - Muller Method.
// input: expectation E and variance V.
// output: random number in Gaussian distribution
double RandG(double E, double V)
{
	static double v, fac;
	static int  phase = 0;

	double S, Z, U1, U2, u;

	if (phase)
		Z = v * fac;
	else
	{
		do
		{
			U1 = RandR();
			U2 = RandR();

			u = 2 * U1 - 1;
			v = 2 * U2 - 1;
			S = u * u + v * v;
		}
		while (S >= 1);

		fac = sqrt(-2 * log(S) / S);
		Z = u * fac;
	}

	phase = 1 - phase;

	return Z;
}
