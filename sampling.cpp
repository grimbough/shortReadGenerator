/*
 * sampling.cpp
 *
 *  Created on: 5 Oct 2011
 *      Author: smith68
 */

#include "shortReadGenerator.h"

int SampleNormal (double mean, double sigma, int seed)
{
    // Create a Mersenne twister random number generator
    static mt19937 rng(static_cast<unsigned> (seed));

    // select Gaussian probability distribution
    normal_distribution<double> norm_dist(mean, sigma);

    // bind random number generator to distribution, forming a function
    variate_generator<mt19937&, normal_distribution<double> >  normal_sampler(rng, norm_dist);

    // sample from the distribution
    return (int) normal_sampler();
}

int SampleUniform (int min, int max, int seed)
{
    static mt19937 rng(static_cast<unsigned> (seed));
    uniform_int<> uni_dist(min, max);
    variate_generator<mt19937&, uniform_int<> >  uniform_sampler(rng, uni_dist);
    return uniform_sampler();
}

double SampleBernoulli (double prob, int seed)
{
    static mt19937 rng(static_cast<unsigned> (seed));
    bernoulli_distribution<> b_dist(prob);
    variate_generator<mt19937&, bernoulli_distribution<double> >  b_sampler(rng, b_dist);
    return b_sampler();
}


