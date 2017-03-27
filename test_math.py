#!/usr/bin/env python3

import numpy as np
import subprocess
import random

import testing




def gemv_test(element_type, tol):
    print("[*] Beginning gemv_test<" + element_type + ">...")
    n = random.randint(5, 1000)
    m = random.randint(5, 1000)

    print("\t n = {}".format(n))
    print("\t m = {}".format(m))

    a = testing.random_matrix(element_type, n, m)
    x = testing.random_matrix(element_type, m, 1)

    y1 = subprocess.check_output(
        ["./tests/gemv", element_type, str(n), str(m)],
        input=(testing.to_cxx_input(a, element_type) 
               + testing.to_cxx_input(x, element_type)
              ).encode('ascii')
    ).decode('ascii').strip('\n')
    y1 = testing.to_matrix(y1, element_type)

    y2 = np.array(np.matrix(a) * np.matrix(x))

    for i in range(n):
        if np.absolute((y1[i] - y2[i]) / y1[i]) > tol:
            raise Exception("Test failed:\n"
                            + "C++ != Python\n"
                            + str(y1[i]) + " != " + str(y2[i])
                           )
    print("[+] Succes!")


def gemm_test(element_type, tol):
    print("[*] Beginning gemm_test<" + element_type + ">...")
    n = random.randint(5, 1000)
    m = random.randint(5, 1000)
    k = random.randint(5, 1000)

    print("\t n = {}".format(n))
    print("\t m = {}".format(m))
    print("\t k = {}".format(k))


    a = testing.random_matrix(element_type, n, m)
    b = testing.random_matrix(element_type, m, k)

    c1 = subprocess.check_output(
        ["./tests/gemm", element_type, str(n), str(m), str(k)],
        input=(testing.to_cxx_input(a, element_type) 
               + testing.to_cxx_input(b, element_type)
              ).encode('ascii') 
    ).decode('ascii').strip('\n')
    c1 = testing.to_matrix(c1, element_type)

    c2 = np.array(np.matrix(a) * np.matrix(b))

    for i in range(n):
        for j in range(k):
            if np.absolute((c1[i, j] - c2[i, j]) / c1[i, j]) > tol:
                raise Exception("Test failed:\n" 
                                + "C++ != Python\n"
                                + str(c1[i, j]) + " != " + str(c2[i, j])
                               )
    print("[+] Succes!")



def dot_test(element_type, tol):
    print("[*] Beginning dot_test<" + element_type + ">...")
    n = random.randint(5, 10000)
    print("\tn = {}".format(n))

    x = testing.random_matrix(element_type, n, 1)
    y = testing.random_matrix(element_type, n, 1)

    r1 = subprocess.check_output(
        ["./tests/dot", element_type, str(n)],
        input=(testing.to_cxx_input(x, element_type) 
               + testing.to_cxx_input(y, element_type)
              ).encode('ascii') 
    ).decode('ascii').strip('\n')

    r1 = testing.to_number(r1, element_type)
    r2 = np.vdot(x, y)

    if np.absolute((r1 - r2) / r1) > tol:
        raise Exception("Test failed!\n"
                        + "C++ != Python\n"
                        + str(r1) + " != " + str(r2)
                       )
    print('[+] Succes!')

def heevr_test(element_type, tol):
    print("[*] Beginning heevr_test<" + element_type + ">...")
    n = 1200 # random.randint(1000, 3000)
    print("\tn = {}".format(n))

    a = testing.random_matrix(element_type, n, n)
    a = a + a.T
    upper_triangle = np.triu_indices(n)
    a[upper_triangle] = np.conjugate(a[upper_triangle])

    w1 = subprocess.check_output( 
        ["./tests/heevr", element_type, str(n)],
        input=(testing.to_cxx_input(a, element_type)
              ).encode('ascii') 
        ).decode('ascii').strip('\n')

    # eigenvalues of a Hermitian matrix are real,
    # thus we remove the complex part of the type
    # if it exists
    if element_type == 'complex-float' \
       or element_type == 'complex-double':
        element_type = element_type[8:]
    w1 = testing.to_matrix(w1, element_type)

    w2 = np.linalg.eigvalsh(a)

    for i in range(n):
        if np.absolute((w1[i] - w2[i]) / w1[i]) > tol:
            raise Exception("Test failed!\n"
                            + "C++ != Python\n"
                            + str(w1[i]) + " != " + str(w2[i])
                           )
    print('Succes!')



def fft_1d_test(element_type, tol):
    print("[*] Beginning fft_1d_test<" + element_type + ">...", end='')

    if element_type == 'float' or element_type == 'double':
        print("Nothing to be done.")
        return

    n = random.randint(10, 1000)
    a = random_matrix(element_type, n, 1)

    b1 = subprocess.check_output(
        ["./fft", element_type, str(n)],
        input=to_cxx_input(a, element_type).encode('ascii')
        ).decode('ascii').strip('\n')

    b1 = to_matrix(b1, element_type)
    b2 = np.fft.fft(a, axis=0) 

    for i in range(n):
        if np.absolute(b1[i, 0] - b2[i, 0]) > tol:
            raise Exception("Test failed!\n" 
                            + str(b1[i]) + " =/= " + str(b2[i])
                           )
    print('Succes!')


            
def fft_2d_test(element_type, tol):
    print("[*] Beginning fft_2d_test<" + element_type + ">...", end='')

    if element_type == 'float' or element_type == 'double':
        print("Nothing to be done.")
        return

    n = random.randint(10, 100)
    m = random.randint(10, 100)

    a = random_matrix(element_type, n, m)

    b1 = subprocess.check_output(
        ["./fft_2d", element_type, str(n), str(m)],
        input=to_cxx_input(a, element_type).encode('ascii')
        ).decode('ascii').strip('\n')
    b1 = to_matrix(b1, element_type)
    
    b2 = np.fft.fft2(a) 

    for i in range(n):
        for j in range(m):
            if np.absolute(b1[i, j] - b2[i, j]) > tol:
                raise Exception("Test failed!\n" 
                                + str(b1[i, j]) + " =/= " + str(b2[i, j])
                               )
    print('Succes!')



def main():

    tests = [heevr_test, 
            # dot_test,
            ]
    types = ['float', 'complex-float', 'double', 'complex-double']
    tol_map = { 'float' : 1.E-3
              , 'double' : 1.E-5
              , 'complex-float' : 1.E-2
              , 'complex-double' : 1.E-5
              }

    for i in range(2):
        print('#' * 30)
        print('# PASS ' + str(i))
        print('#' * 30)
        for test in tests:
            for element_type in types:
                test(element_type, tol_map[element_type])


if __name__ == '__main__':
    main()
