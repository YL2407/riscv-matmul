.equ r1, 4
.equ c1, 2
#c1 = r2 by assumption
.equ c2, 3

#arr 1:
#40
#-10
#1.235
#-0.734
#2.885
#-7.671
#3.521
#8.12
#arr 2:
#-7.5
#0.02345
#-2.397
#-1.275
#0.1257
#1.591
#-7.73
#1.023
.data
first_array:
    .word 0x42200000, 0xc1200000
    .word 0x3f9e147b, 0xbf3be76d, 0x4038a3d7, 0xc0f578d5, 0x40615810, 0x4101eb85
second_array:
    #.word 0xC0F00000, 0x3cc01a37
    .word 0xc0196873, 0xbfa33333, 0x3e00b780, 0x3fcba5e3, 0xc0f75c29, 0x3f82f1aa
result_array:
    .word 0
.text
la t0 first_array
la t1 second_array
la tp result_array

#now for the matrix multiplication loop code
mm:
    #we loop over rows of first and columns of second, sum the products (dot product)
    #then we move onto the next
    li s7, r1
    li s8, c1
    li s9, c2
    outer:
        li s9, c2
        la t1 second_array
        beq s7, zero, end
        inner:
            beq s9, zero, next_outer
            li s10, 0
            li s11, 0
            li a0, 0
            dp:
            beq s10, s8, next_inner #don't forget to implement this
            slli s10, s10, 2#remember to undo
            add s10, s10, t0
            slli s11, s11, 2
            add s11, s11, t1
            lw s0, 0(s10)
            lw s1, 0(s11)
            j fp_mult
            return:
                #floating point addition into a0
                beq a0, zero, init

                srli a1, a0, 23
                andi a1, a1, 0xFF
                srli a2, s5, 23
                andi a2, a2, 0xFF
                li t6 0x7FFFFF
                li t4, 1
                slli t4, t4, 23
                #exponents now stored in a1 and a2
                blt a1, a2, shift_first
                shift_second:
                    mv a3, a1
                    sub t5, a1, a2 #absolute diff
                    and a1, a0, t6
                    and a2, s5, t6
                    or a1, a1, t4
                    or a2, a2, t4
                    
                    srl a2, a2, t5
                    j add_significands
                shift_first:
                    mv a3, a2
                    sub t5, a2, a1 #absolute diff
                    and a1, a0, t6
                    and a2, s5, t6
                    or a1, a1, t4
                    or a2, a2, t4
                    srl a1, a1, t5
                #a1 and a2 now have significands, appropriately shifted, implicit 1 added back in
                #a3 contains the correct exponent for the unnormalised sum
                add_significands:
                    #ah, based on sign, remember
                    #extract the sign
                    srli t6, a0, 31
                    beq t6, zero, skip_first
                    sub a1, zero, a1
                    skip_first:
                    srli t5, s5, 31
                    beq t5, zero, skip_second
                    sub a2, zero, a2
                    skip_second:
                    add a4, a1, a2
                    blt a4, zero, adjust_sign
                    li a5, 0
                    j normalise_sum
                    adjust_sign:
                    li a5, 1
                    sub a4, zero, a4
                    #a4 contains the unnormalised sum, a5 contains the sign
                normalise_sum:
                    li t6, 8
                    li t5, 1
                    slli t5, t5, 31
                    nloop:
                    and t4, a4, t5
                    bne t4, zero, found_msb
                    addi t6, t6, -1
                    srli t5, t5, 1
                    #again, assuming there is at least one 1
                    j nloop
                found_msb:
                    add a3, a3, t6
                    #now normalised, we just need to shift the significand appropriately
                    #ah, we cant shift by negatives, can we
                    blt t6, zero, shift_left
                    srl a4, a4, t6
                    li t5, 0x7FFFFF
                    and a4, a4, t5 #remove implicit 1
                    j reconstruct
                    shift_left:
                        sub t6, zero, t6
                        sll a4, a4, t6
                        li t5, 0x7FFFFF
                        and a4, a4, t5 #remove implicit 1
                reconstruct:
                    mv a6, a5
                    slli a6, a6, 8
                    andi a3, a3, 0xFF 
                    or a6, a6, a3
                    slli a6, a6, 23
                    or a6, a6, a4
                mv a0, a6
                j continue
                #now a0 contains the sum of a0 and s5 as floating point numbers
            init:
                mv a0, s5
            continue:
                sub s10, s10, t0
                srli s10, s10, 2
                addi s10, s10, 1
                sub s11, s11, t1
                srli s11, s11, 2
                addi s11, s11, c2 
                j dp
            next_inner:
                addi s9, s9, -1
                addi t1, t1, 4
                sw a0, 0(tp)
                addi tp, tp, 4
                j inner
        next_outer:
            addi s7, s7, -1
            li a7, c1
            slli a7, a7, 2 #cols * 4
            add t0, t0, a7
            j outer
end:
    nop
    j absolute_end


#temporary
#lw s0, 0(t0)
#lw s1, 0(t1)

fp_mult: #assuming numbers are stored in s0 and s1
    li s2, 0
    li t4, 1
    slli t4, t4, 23
    li t3, 0x7FFFFF
    and t6, s1, t3
    and t5, s0, t3
    or t6, t6, t4
    or t5, t5, t4
    #now we have extracted the significands and added back in the implicit 1
    li t2, 1
    li s6, 17
    #keeping track of what bit we're at in the multiplication in t2
    loop1:
        beq t6, zero, exponents
        andi t4, t6, 1
        beq t4, zero, skip
        add s2, s2, t5
        skip:
            blt t2, s6, shift_sum
            shift_num:
                slli t5, t5, 1
                j uncond
            shift_sum:
                srli s2, s2, 1
            uncond:
                srli t6, t6, 1
                addi t2, t2, 1
                j loop1
    #s2 now contains the product of the significands
    exponents:
        srli t5, s0, 23
        andi t5, t5, 0xFF
        srli t6, s1, 23
        andi t6, t6, 0xFF
        add t5, t5, t6
        li t6, 1
        slli t6, t6, 7
        addi t6, t6, -1 #bias
        sub t5, t5, t6
        mv s3, t5
    #s3 now contains the sum of the exponents represented correctly
    sign:
        xor t5, s0, s1
        srli t5, t5, 31
        mv s4, t5
    #s4 now contains the sign of the result
    normalise:
        #we want to find the first 1 in the significand and adjust the exponent so the radix point is just after
        #radix point starts off after 2 digits
        #we make the shocking assumption there is at least one 1, i.e., we never multiply by 0 (and indeed, even the fp representation of 0 has a 1)
        li t5, 1 #exponent shift amount
        mv t6, s2
        li t4, 1
        slli t4, t4, 31 #mask to check if MSB is 1
        li t2, 1#getting rid of the implicit 1
        loop2:
            and t3, t6, t4
            bnez t3, found_one
            addi t5, t5, -1
            addi t2, t2, 1
            slli t6, t6, 1
            j loop2
        found_one:
            add s3, s3, t5
            sll s2, s2, t2 #remove implicit 1
        #s3 now contains correct exponent for normalised significand
    combine:
        slli s5, s4, 31
        slli t6, s3, 23
        or s5, s5, t6
        li t5, 1
        slli t5, t5, 23
        addi t5, t5 -1
        #srli s2, s2, 9
        #slli s2, s2 1
        srli s2, s2, 9
        and t5, t5, s2 #remove the implicit 1/mask to fit within the bits (just in case)
        or s5, s5, t5
    j return
    nop

absolute_end:
    nop
