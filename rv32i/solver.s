# solver.s: RV32I IDA* solver for the 2x2x2 cube.
# Same algorithm and move order as search() in solver_new.c.
# An optimal path need not match tests/solutions.txt: several shortest paths can exist.
# Ripes has no .include; build.ps1 joins this file with tables.s.
# All pseudoinstructions here expand to base RV32I instructions.
# This is the CLI search build; the LED renderer is a separate next step.

        .data
input:  .string "21345671111111"    # state to solve: 7 cubies, then 7 twists
        .align 2
expected_length: .word -1           # -1: arbitrary query; >= 0: test assertion
face_names: .string "RBD"           # face 0 = R, face 1 = B, face 2 = D
failure_message: .string "FAIL\n"

        .bss
        .align 2
# One 24-byte record per depth g = 0..11:
#   0: P   permutation rank        4: O   orientation rank
#   8: F   face being tried       12: T   turn being tried
#  16: PB  byte offset of face F in perm_move (F * 10080)
#  20: OB  byte offset of face F in ori_move  (F * 1458)
frames: .zero 288
lehmer: .zero 7                    # scratch space for the seven Lehmer digits

# Register map (no calls, no recursion, no machine stack needed):
# s0, s1: original P, O, preserved until the independent T5 replay
# s2, s3: perm_move, ori_move base addresses
# s4, s5: perm_dist, ori_dist base addresses
# s6: current frame address = frames + 24*g
# s7: current depth g; also the actual solution length after success
# s8: IDA* bound; s9: 3; s10: 11; s11: permutation byte stride 10080
# t0..t6: scratch registers; a0/a7: Ripes syscall arguments

        .text
main:
# ---- 1. Input -> two ranks (parse_ranks in C) -------------------------
# ASCII comparison gives the same ordering as subtracting '1' first.
        la t0, input               # current cubie character
        addi t1, t0, 7             # end of the cubie characters
        la t2, lehmer              # destination for c[i]
parse_cubie:
        lbu t3, 0(t0)
        addi t4, t0, 1             # first later character, j = i+1
        li t5, 0                   # smaller = 0
count_smaller:
        beq t4, t1, store_digit
        lbu t6, 0(t4)
        sltu t6, t6, t3            # 1 if the later cubie is smaller
        add t5, t5, t6
        addi t4, t4, 1
        j count_smaller
store_digit:
        sb t5, 0(t2)
        addi t2, t2, 1
        addi t0, t0, 1
        bne t0, t1, parse_cubie

# Horner: (((((c[0]*6+c[1])*5+c[2])*4+c[3])*3+c[4])*2+c[5]).
# c[6] is always zero, so the final '*1 + c[6]' changes nothing.
        la t0, lehmer
        lbu s0, 0(t0)
        slli t1, s0, 2
        slli t2, s0, 1
        add s0, t1, t2             # r*6 = (r<<2) + (r<<1)
        lbu t3, 1(t0)
        add s0, s0, t3
        slli t1, s0, 2
        add s0, t1, s0             # r*5 = (r<<2) + r
        lbu t3, 2(t0)
        add s0, s0, t3
        slli s0, s0, 2             # r*4
        lbu t3, 3(t0)
        add s0, s0, t3
        slli t1, s0, 1
        add s0, t1, s0             # r*3 = (r<<1) + r
        lbu t3, 4(t0)
        add s0, s0, t3
        slli s0, s0, 1             # r*2
        lbu t3, 5(t0)
        add s0, s0, t3

# The seventh twist is implied by sum(twists) == 0 (mod 3) for valid input.
        la t0, input
        addi t0, t0, 7
        li t1, 6
        li s1, 0
parse_orientation:
        lbu t2, 0(t0)
        addi t2, t2, -49           # digit - '1'
        slli t3, s1, 1
        add s1, t3, s1             # q*3
        add s1, s1, t2
        addi t0, t0, 1
        addi t1, t1, -1
        bnez t1, parse_orientation

# ---- 2. Iterative deepening (solve in C) ------------------------------
        la s2, perm_move
        la s3, ori_move
        la s4, perm_dist
        la s5, ori_dist
        li s9, 3
        li s10, 11
        li s11, 10080              # 5040 halfwords, in bytes
        add t0, s4, s0
        lbu s8, 0(t0)
        add t0, s5, s1
        lbu t1, 0(t0)
        bge s8, t1, iterate
        mv s8, t1                 # bound = max(perm_dist[P], ori_dist[O])
iterate:
        la s6, frames
        li s7, 0
        sw s0, 0(s6)
        sw s1, 4(s6)

# ---- 3. search: one IDA* iteration (search() in C) ---------------------
enter:
        lw t0, 0(s6)               # P[g]
        lw t1, 4(s6)               # O[g]
        or t2, t0, t1
        beqz t2, print             # P == 0 and O == 0
# g + max(hp, ho) > bound iff hp > bound-g or ho > bound-g.
# Test permutation first: a failure avoids the orientation table load.
        sub t2, s8, s7             # moves remaining in this iteration
        add t3, s4, t0
        lbu t3, 0(t3)
        blt t2, t3, back
        add t3, s5, t1
        lbu t3, 0(t3)
        blt t2, t3, back
# Only (0,0) has heuristic zero, so an unsolved node at g==bound is
# always pruned above; no child can write past frame 11.
        li t2, -1
        sw t2, 8(s6)               # F = -1
        sub t2, zero, s11
        sw t2, 16(s6)              # PB = -10080
        li t2, -1458
        sw t2, 20(s6)              # OB = -1458
next_face:
        lw t0, 8(s6)
        addi t0, t0, 1
        sw t0, 8(s6)               # F++
        beq t0, s9, back           # all three faces have been tried
        lw t1, 16(s6)
        add t1, t1, s11
        sw t1, 16(s6)              # PB += 10080
        lw t2, 20(s6)
        addi t2, t2, 1458
        sw t2, 20(s6)              # OB += 1458
        beqz s7, first_turn        # root has no parent
        lw t1, -16(s6)             # parent F: -24 + 8
        beq t0, t1, next_face      # consecutive moves of one face combine
first_turn:
        sw zero, 12(s6)            # T = 0 (one quarter turn)
        lw t0, 0(s6)
        lw t1, 4(s6)
        j turn_lookup
next_turn:
        lw t0, 12(s6)
        addi t0, t0, 1
        beq t0, s9, next_face      # T == 3: this face is done
        sw t0, 12(s6)
        lw t0, 24(s6)              # previous child P, NOT parent P
        lw t1, 28(s6)              # previous child O
turn_lookup:
# perm_move/ori_move are uint16_t: the rank must be shifted by one.
# PB/OB already contain BYTE offsets, unlike the C element offsets.
        lw t2, 16(s6)
        add t2, s2, t2
        slli t0, t0, 1
        add t0, t2, t0
        lhu t0, 0(t0)
        sw t0, 24(s6)              # child P
        lw t2, 20(s6)
        add t2, s3, t2
        slli t1, t1, 1
        add t1, t2, t1
        lhu t1, 0(t1)
        sw t1, 28(s6)              # child O
child:
        addi s7, s7, 1
        addi s6, s6, 24
        j enter
back:
        beqz s7, iteration_failed
        addi s7, s7, -1
        addi s6, s6, -24
        j next_turn
iteration_failed:
        addi s8, s8, 1
        blt s10, s8, fail          # bound > 11
        j iterate

# ---- 4. Print the moves ------------------------------------------------
print:
# The successful frames are the path: no separate path[] array needed.
# Use the actual depth s7 for both printing and replay.
        la t0, expected_length
        lw t0, 0(t0)
        blt t0, zero, print_begin
        bne t0, s7, fail           # optional in-program optimal-length check
print_begin:
        la t0, frames
        li t1, 0
        la t3, face_names
        li a7, 11                 # Ripes print-character syscall
print_move:
        beq t1, s7, replay_begin
        lw t2, 8(t0)
        add t2, t3, t2
        lbu a0, 0(t2)
        ecall
        lw t2, 12(t0)
        beqz t2, print_space       # T == 0: R/B/D
        li t4, 1
        li a0, 50                 # '2' for T == 1
        beq t2, t4, print_suffix
        li a0, 39                 # apostrophe for T == 2
print_suffix:
        ecall
print_space:
        li a0, 32
        ecall
        addi t0, t0, 24
        addi t1, t1, 1
        j print_move

# ---- 5. T5: apply the moves to (p, o) and check (0, 0) ----------------
replay_begin:
        la t0, frames
        li t1, 0
replay_move:
        beq t1, s7, replay_done
        lw t2, 16(t0)
        add t2, s2, t2             # selected face's permutation table base
        lw t3, 20(t0)
        add t3, s3, t3             # selected face's orientation table base
        lw t4, 12(t0)
        addi t4, t4, 1             # T=0,1,2 means 1,2,3 quarter turns
replay_quarter:
        slli t5, s0, 1
        add t5, t2, t5
        lhu s0, 0(t5)
        slli t5, s1, 1
        add t5, t3, t5
        lhu s1, 0(t5)
        addi t4, t4, -1
        bnez t4, replay_quarter
        addi t0, t0, 24
        addi t1, t1, 1
        j replay_move
replay_done:
        li a0, 10                 # newline (also for the solved input)
        li a7, 11
        ecall
        or a0, s0, s1
        snez a0, a0               # exit status 0 iff replay reached (0,0)
        li a7, 93
        ecall
fail:
        la a0, failure_message
        li a7, 4
        ecall
        li a0, 1
        li a7, 93
        ecall
