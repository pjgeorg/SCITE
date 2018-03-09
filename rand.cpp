/*
 * rand.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 *  Modified on: Mar, 2018
 *      Author: pjgeorg
 */

#include "rand.h"

std::default_random_engine RandomGenerator::sEngine((std::random_device())());
