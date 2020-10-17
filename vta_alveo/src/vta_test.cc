/*!
 *  Copyright (c) 2018 by Contributors
 * \file accel_test.cpp
 * \brief Simulation tests for the VTA design.
 */

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "./test_lib.h"

int main(void) {

#if VTA_DEBUG == 1
    printParameters();
#endif

    // Run load/store test
    int status = mem_test(128, 128);
    // Run reset test
    status |= reset_test(128, 128);
    // Run fully connected layer test
    status |= fc_test(64, 128, 128);

    // Run 2D convolution layer test
    // TODO Part 2: Uncomment line below
    // status |= conv2d_test(1, 9, 9, 3, 3, 32, 64);

    return status;
}
