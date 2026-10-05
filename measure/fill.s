# Measurement 1
# Get to know the used time of hash map in Ripes
    .equ BASE, 0x10000000 #starting point
    .equ COUNT, 65536
    .equ STRIDE, 4

    .text
main:    
    li t0, BASE
    li t1, COUNT
    li t2, STRIDE
loop:
    sw zero, 0(t0)
    add t0, t0, t2
    addi t1, t1, -1
    bnez t1, loop #go back to loop until t1 is zero
    li a7, 10
    ecall