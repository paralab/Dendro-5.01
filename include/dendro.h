//
// Created by milinda on 2/8/16.
//

#ifndef SFCSORTBENCH_DENDRO_H
#define SFCSORTBENCH_DENDRO_H

#define RED "\e[1;31m"
#define BLU "\e[2;34m"
#define GRN "\e[0;32m"
#define YLW "\e[0;33m"
#define MAG "\e[0;35m"
#define CYN "\e[0;36m"
#define NRM "\e[0m"



#ifdef USE_64BIT_INDICES
#define DendroIntL long long
#define DendroIntLSpecifier %lld
#define DendroUIntLSpecifier %llu
#else
#define DendroIntL unsigned int
#define DendroIntLSpecifier %d
#define DendroUIntLSpecifier %u
#endif


//#define DendroIntL unsigned int
typedef unsigned __int128 DendroUInt_128;

#endif //SFCSORTBENCH_DENDRO_H
