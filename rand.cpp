/*
 * rand_C.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

//#include <string>
//#include <random>
#include "rand.h"
#include <iostream>
#include <random>
#include <stdlib.h>
#include <time.h>
#include "matrices.h"

std::default_random_engine RandomGenerator::sEngine((std::random_device())());
