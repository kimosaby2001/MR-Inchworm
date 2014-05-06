
;; Function kmer_int_type_t revcomp_val(kmer_int_type_t, unsigned int) (_Z11revcomp_valyj)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
merging block 8 into block 7
Merged 7 and 8 without moving.
merging block 9 into block 7
Merged 7 and 9 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 3 0 7 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 7 3 4 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 4 7 5 2 sequenceUtil.cpp:182 (set (reg/v:DI 65 [ kmer ])
        (reg:DI 5 di [ kmer ])) -1 (nil))

(insn 5 4 6 2 sequenceUtil.cpp:182 (set (reg/v:SI 66 [ kmer_length ])
        (reg:SI 4 si [ kmer_length ])) -1 (nil))

(note 6 5 9 2 NOTE_INSN_FUNCTION_BEG)

(debug_insn 9 6 10 2 sequenceUtil.cpp:183 (var_location:DI rev_kmer (const_int 0 [0x0])) -1 (nil))

(debug_insn 10 9 11 2 sequenceUtil.cpp:184 (var_location:DI D.4294967180 (not:DI (reg/v:DI 65 [ kmer ]))) -1 (nil))

(debug_insn 11 10 12 2 sequenceUtil.cpp:184 (var_location:DI kmer (debug_expr:DI D#116)) -1 (nil))

(debug_insn 12 11 13 2 sequenceUtil.cpp:185 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 13 12 14 2 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 14 13 15 2 (var_location:DI rev_kmer (const_int 0 [0x0])) -1 (nil))

(debug_insn 15 14 16 2 (var_location:DI kmer (debug_expr:DI D#116)) -1 (nil))

(insn 16 15 17 2 sequenceUtil.cpp:185 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 66 [ kmer_length ])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 17 16 18 2 sequenceUtil.cpp:185 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 22)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9100 [0x238c])
        (nil)))
;; End of basic block 2 -> ( 4 3)

;; Succ edge  4 [91.0%] 
;; Succ edge  3 [9.0%]  (fallthru)

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [9.0%]  (fallthru)
(note 18 17 19 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 19 18 20 3 sequenceUtil.cpp:185 (set (reg/v:DI 60 [ rev_kmer.537 ])
        (const_int 0 [0x0])) -1 (nil))

(jump_insn 20 19 21 3 sequenceUtil.cpp:185 (set (pc)
        (label_ref 45)) -1 (nil))
;; End of basic block 3 -> ( 6)

;; Succ edge  6 [100.0%] 

(barrier 21 20 22)

;; Start of basic block ( 2) -> 4
;; Pred edge  2 [91.0%] 
(code_label 22 21 23 4 2 "" [1 uses])

(note 23 22 24 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 24 23 25 4 sequenceUtil.cpp:184 (set (reg/v:DI 61 [ kmer.536 ])
        (not:DI (reg/v:DI 65 [ kmer ]))) -1 (nil))

(insn 25 24 26 4 sequenceUtil.cpp:184 (set (reg/v:SI 62 [ i ])
        (const_int 0 [0x0])) -1 (nil))

(insn 26 25 42 4 sequenceUtil.cpp:184 (set (reg/v:DI 60 [ rev_kmer.537 ])
        (const_int 0 [0x0])) -1 (nil))
;; End of basic block 4 -> ( 5)

;; Succ edge  5 [100.0%]  (fallthru)

;; Start of basic block ( 5 4) -> 5
;; Pred edge  5 [91.0%]  (dfs_back)
;; Pred edge  4 [100.0%]  (fallthru)
(code_label 42 26 27 5 4 "" [1 uses])

(note 27 42 28 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(debug_insn 28 27 29 5 sequenceUtil.cpp:187 (var_location:SI base (and:SI (subreg:SI (reg/v:DI 61 [ kmer.536 ]) 0)
        (const_int 3 [0x3]))) -1 (nil))

(insn 29 28 30 5 sequenceUtil.cpp:188 (parallel [
            (set (reg/v:DI 63 [ rev_kmer ])
                (ashift:DI (reg/v:DI 60 [ rev_kmer.537 ])
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 30 29 31 5 sequenceUtil.cpp:188 (var_location:DI rev_kmer (reg/v:DI 63 [ rev_kmer ])) -1 (nil))

(insn 31 30 32 5 sequenceUtil.cpp:189 (parallel [
            (set (reg:SI 67)
                (and:SI (subreg:SI (reg/v:DI 61 [ kmer.536 ]) 0)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 32 31 33 5 sequenceUtil.cpp:189 (set (reg:DI 68)
        (sign_extend:DI (reg:SI 67))) -1 (nil))

(insn 33 32 34 5 sequenceUtil.cpp:189 (parallel [
            (set (reg/v:DI 60 [ rev_kmer.537 ])
                (plus:DI (reg:DI 68)
                    (reg/v:DI 63 [ rev_kmer ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 34 33 35 5 sequenceUtil.cpp:189 (var_location:DI rev_kmer (reg/v:DI 60 [ rev_kmer.537 ])) -1 (nil))

(insn 35 34 36 5 sequenceUtil.cpp:190 (parallel [
            (set (reg/v:DI 61 [ kmer.536 ])
                (lshiftrt:DI (reg/v:DI 61 [ kmer.536 ])
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 36 35 37 5 sequenceUtil.cpp:190 (var_location:DI kmer (reg/v:DI 61 [ kmer.536 ])) -1 (nil))

(insn 37 36 38 5 sequenceUtil.cpp:185 (parallel [
            (set (reg/v:SI 62 [ i ])
                (plus:SI (reg/v:SI 62 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 38 37 39 5 sequenceUtil.cpp:185 (var_location:SI i (reg/v:SI 62 [ i ])) -1 (nil))

(debug_insn 39 38 40 5 (var_location:SI i (reg/v:SI 62 [ i ])) -1 (nil))

(debug_insn 40 39 41 5 (var_location:DI rev_kmer (reg/v:DI 60 [ rev_kmer.537 ])) -1 (nil))

(debug_insn 41 40 43 5 (var_location:DI kmer (reg/v:DI 61 [ kmer.536 ])) -1 (nil))

(insn 43 41 44 5 sequenceUtil.cpp:185 (set (reg:CC 17 flags)
        (compare:CC (reg/v:SI 66 [ kmer_length ])
            (reg/v:SI 62 [ i ]))) -1 (nil))

(jump_insn 44 43 45 5 sequenceUtil.cpp:185 (set (pc)
        (if_then_else (gtu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 42)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9100 [0x238c])
        (nil)))
;; End of basic block 5 -> ( 5 6)

;; Succ edge  5 [91.0%]  (dfs_back)
;; Succ edge  6 [9.0%]  (fallthru)

;; Start of basic block ( 5 3) -> 6
;; Pred edge  5 [9.0%]  (fallthru)
;; Pred edge  3 [100.0%] 
(code_label 45 44 46 6 3 "" [1 uses])

(note 46 45 47 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 47 46 51 6 sequenceUtil.cpp:185 (set (reg:DI 64 [ <result> ])
        (reg/v:DI 60 [ rev_kmer.537 ])) -1 (nil))

(insn 51 47 57 6 sequenceUtil.cpp:195 (set (reg/i:DI 0 ax)
        (reg:DI 64 [ <result> ])) -1 (nil))

(insn 57 51 0 6 sequenceUtil.cpp:195 (use (reg/i:DI 0 ax)) -1 (nil))
;; End of basic block 6 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function int base_to_int_value(char) (_Z17base_to_int_valuec)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
merging block 7 into block 6
Merged 6 and 7 without moving.
merging block 8 into block 6
Merged 6 and 8 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 5 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 5 1 2 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 2 5 3 2 sequenceUtil.cpp:213 (set (reg:SI 62)
        (reg:SI 5 di [ nucleotide ])) -1 (nil))

(insn 3 2 4 2 sequenceUtil.cpp:213 (set (reg/v:QI 61 [ nucleotide ])
        (subreg:QI (reg:SI 62) 0)) -1 (nil))

(note 4 3 7 2 NOTE_INSN_FUNCTION_BEG)

(insn 7 4 8 2 sequenceUtil.cpp:213 (parallel [
            (set (reg:QI 58 [ csui.434 ])
                (plus:QI (reg/v:QI 61 [ nucleotide ])
                    (const_int -65 [0xffffffffffffffbf])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 8 7 9 2 sequenceUtil.cpp:213 (set (reg:CC 17 flags)
        (compare:CC (reg:QI 58 [ csui.434 ])
            (const_int 51 [0x33]))) -1 (nil))

(jump_insn 9 8 10 2 sequenceUtil.cpp:213 (set (pc)
        (if_then_else (leu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 14)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 6100 [0x17d4])
        (nil)))
;; End of basic block 2 -> ( 3 4)

;; Succ edge  3 [39.0%]  (fallthru)
;; Succ edge  4 [61.0%] 

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [39.0%]  (fallthru)
(note 10 9 11 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 11 10 12 3 sequenceUtil.cpp:213 (set (reg:SI 59 [ D.32203 ])
        (const_int -1 [0xffffffffffffffff])) -1 (nil))

(jump_insn 12 11 13 3 sequenceUtil.cpp:213 (set (pc)
        (label_ref 19)) -1 (nil))
;; End of basic block 3 -> ( 5)

;; Succ edge  5 [100.0%] 

(barrier 13 12 14)

;; Start of basic block ( 2) -> 4
;; Pred edge  2 [61.0%] 
(code_label 14 13 15 4 9 "" [1 uses])

(note 15 14 16 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 16 15 17 4 sequenceUtil.cpp:213 (set (reg:DI 63)
        (zero_extend:DI (reg:QI 58 [ csui.434 ]))) -1 (nil))

(insn 17 16 18 4 sequenceUtil.cpp:213 (set (reg/f:DI 64)
        (symbol_ref:DI ("CSWTCH.435") [flags 0x2]  <var_decl 0x7f9011cde280 CSWTCH.435>)) -1 (nil))

(insn 18 17 19 4 sequenceUtil.cpp:213 (set (reg:SI 59 [ D.32203 ])
        (mem/s/u:SI (plus:DI (mult:DI (reg:DI 63)
                    (const_int 4 [0x4]))
                (reg/f:DI 64)) [5 CSWTCH.435 S4 A32])) -1 (nil))
;; End of basic block 4 -> ( 5)

;; Succ edge  5 [100.0%]  (fallthru)

;; Start of basic block ( 3 4) -> 5
;; Pred edge  3 [100.0%] 
;; Pred edge  4 [100.0%]  (fallthru)
(code_label 19 18 20 5 10 "" [1 uses])

(note 20 19 21 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 21 20 25 5 sequenceUtil.cpp:213 (set (reg:SI 60 [ <result> ])
        (reg:SI 59 [ D.32203 ])) -1 (nil))

(insn 25 21 31 5 sequenceUtil.cpp:238 (set (reg/i:SI 0 ax)
        (reg:SI 60 [ <result> ])) -1 (nil))

(insn 31 25 0 5 sequenceUtil.cpp:238 (use (reg/i:SI 0 ax)) -1 (nil))
;; End of basic block 5 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function kmer_int_type_t get_maximum_kmer_intval(unsigned int) (_Z23get_maximum_kmer_intvalj)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
merging block 8 into block 7
Merged 7 and 8 without moving.
merging block 9 into block 7
Merged 7 and 9 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 2 0 5 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 5 2 3 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 3 5 4 2 sequenceUtil.cpp:250 (set (reg/v:SI 64 [ kmer_length ])
        (reg:SI 5 di [ kmer_length ])) -1 (nil))

(note 4 3 7 2 NOTE_INSN_FUNCTION_BEG)

(debug_insn 7 4 8 2 sequenceUtil.cpp:252 (var_location:QI c (const_int 67 [0x43])) -1 (nil))

(debug_insn 8 7 9 2 sequenceUtil.cpp:253 (var_location:DI max_kmer (const_int 0 [0x0])) -1 (nil))

(debug_insn 9 8 10 2 sequenceUtil.cpp:254 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 10 9 11 2 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 11 10 12 2 (var_location:DI max_kmer (const_int 0 [0x0])) -1 (nil))

(insn 12 11 13 2 sequenceUtil.cpp:254 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 64 [ kmer_length ])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 13 12 14 2 sequenceUtil.cpp:254 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 18)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9100 [0x238c])
        (nil)))
;; End of basic block 2 -> ( 4 3)

;; Succ edge  4 [91.0%] 
;; Succ edge  3 [9.0%]  (fallthru)

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [9.0%]  (fallthru)
(note 14 13 15 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 15 14 16 3 sequenceUtil.cpp:254 (set (reg/v:DI 59 [ max_kmer.591 ])
        (const_int 0 [0x0])) -1 (nil))

(jump_insn 16 15 17 3 sequenceUtil.cpp:254 (set (pc)
        (label_ref 38)) -1 (nil))
;; End of basic block 3 -> ( 6)

;; Succ edge  6 [100.0%] 

(barrier 17 16 18)

;; Start of basic block ( 2) -> 4
;; Pred edge  2 [91.0%] 
(code_label 18 17 19 4 13 "" [1 uses])

(note 19 18 20 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 20 19 21 4 sequenceUtil.cpp:254 (set (reg/f:DI 65)
        (const:DI (plus:DI (symbol_ref:DI ("_base_to_int") [flags 0x2]  <var_decl 0x7f9012470be0 _base_to_int>)
                (const_int 67 [0x43])))) -1 (nil))

(insn 21 20 22 4 sequenceUtil.cpp:254 (parallel [
            (set (reg:SI 66)
                (zero_extend:SI (mem/s/j:QI (reg/f:DI 65) [0 _base_to_int+67 S1 A8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 22 21 23 4 sequenceUtil.cpp:254 (set (reg:DI 60 [ pretmp.576 ])
        (sign_extend:DI (reg:SI 66))) -1 (nil))

(insn 23 22 24 4 sequenceUtil.cpp:254 (set (reg/v:SI 61 [ i ])
        (const_int 0 [0x0])) -1 (nil))

(insn 24 23 35 4 sequenceUtil.cpp:254 (set (reg/v:DI 59 [ max_kmer.591 ])
        (const_int 0 [0x0])) -1 (nil))
;; End of basic block 4 -> ( 5)

;; Succ edge  5 [100.0%]  (fallthru)

;; Start of basic block ( 5 4) -> 5
;; Pred edge  5 [91.0%]  (dfs_back)
;; Pred edge  4 [100.0%]  (fallthru)
(code_label 35 24 25 5 15 "" [1 uses])

(note 25 35 26 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(debug_insn 26 25 27 5 sequenceUtil.cpp:255 (var_location:SI val (zero_extend:SI (mem/s/j:QI (const:DI (plus:DI (symbol_ref:DI ("_base_to_int") [flags 0x2]  <var_decl 0x7f9012470be0 _base_to_int>)
                    (const_int 67 [0x43]))) [0 _base_to_int+67 S1 A8]))) -1 (nil))

(insn 27 26 28 5 sequenceUtil.cpp:256 (parallel [
            (set (reg/v:DI 62 [ max_kmer ])
                (ashift:DI (reg/v:DI 59 [ max_kmer.591 ])
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 28 27 29 5 sequenceUtil.cpp:256 (var_location:DI max_kmer (reg/v:DI 62 [ max_kmer ])) -1 (nil))

(insn 29 28 30 5 sequenceUtil.cpp:257 (parallel [
            (set (reg/v:DI 59 [ max_kmer.591 ])
                (ior:DI (reg/v:DI 62 [ max_kmer ])
                    (reg:DI 60 [ pretmp.576 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 30 29 31 5 sequenceUtil.cpp:257 (var_location:DI max_kmer (reg/v:DI 59 [ max_kmer.591 ])) -1 (nil))

(insn 31 30 32 5 sequenceUtil.cpp:254 (parallel [
            (set (reg/v:SI 61 [ i ])
                (plus:SI (reg/v:SI 61 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 32 31 33 5 sequenceUtil.cpp:254 (var_location:SI i (reg/v:SI 61 [ i ])) -1 (nil))

(debug_insn 33 32 34 5 (var_location:SI i (reg/v:SI 61 [ i ])) -1 (nil))

(debug_insn 34 33 36 5 (var_location:DI max_kmer (reg/v:DI 59 [ max_kmer.591 ])) -1 (nil))

(insn 36 34 37 5 sequenceUtil.cpp:254 (set (reg:CC 17 flags)
        (compare:CC (reg/v:SI 64 [ kmer_length ])
            (reg/v:SI 61 [ i ]))) -1 (nil))

(jump_insn 37 36 38 5 sequenceUtil.cpp:254 (set (pc)
        (if_then_else (gtu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 35)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9100 [0x238c])
        (nil)))
;; End of basic block 5 -> ( 5 6)

;; Succ edge  5 [91.0%]  (dfs_back)
;; Succ edge  6 [9.0%]  (fallthru)

;; Start of basic block ( 5 3) -> 6
;; Pred edge  5 [9.0%]  (fallthru)
;; Pred edge  3 [100.0%] 
(code_label 38 37 39 6 14 "" [1 uses])

(note 39 38 40 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 40 39 44 6 sequenceUtil.cpp:254 (set (reg:DI 63 [ <result> ])
        (reg/v:DI 59 [ max_kmer.591 ])) -1 (nil))

(insn 44 40 50 6 sequenceUtil.cpp:260 (set (reg/i:DI 0 ax)
        (reg:DI 63 [ <result> ])) -1 (nil))

(insn 50 44 0 6 sequenceUtil.cpp:260 (use (reg/i:DI 0 ax)) -1 (nil))
;; End of basic block 6 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function kmer_int_type_t get_DS_kmer_val(kmer_int_type_t, unsigned int) (_Z15get_DS_kmer_valyj)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
merging block 8 into block 7
Merged 7 and 8 without moving.
merging block 9 into block 7
Merged 7 and 9 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 2 0 6 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 6 2 3 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 2 sequenceUtil.cpp:396 (set (reg/v:DI 64 [ kmer_val ])
        (reg:DI 5 di [ kmer_val ])) -1 (nil))

(insn 4 3 5 2 sequenceUtil.cpp:396 (set (reg/v:SI 65 [ kmer_length ])
        (reg:SI 4 si [ kmer_length ])) -1 (nil))

(note 5 4 8 2 NOTE_INSN_FUNCTION_BEG)

(debug_insn 8 5 9 2 (var_location:DI kmer (reg/v:DI 64 [ kmer_val ])) -1 (nil))

(debug_insn 9 8 10 2 (var_location:SI kmer_length (reg/v:SI 65 [ kmer_length ])) -1 (nil))

(debug_insn 10 9 11 2 sequenceUtil.cpp:183 (var_location:DI rev_kmer (const_int 0 [0x0])) -1 (nil))

(debug_insn 11 10 12 2 sequenceUtil.cpp:184 (var_location:DI D.4294967179 (not:DI (reg/v:DI 64 [ kmer_val ]))) -1 (nil))

(debug_insn 12 11 13 2 sequenceUtil.cpp:184 (var_location:DI kmer (debug_expr:DI D#117)) -1 (nil))

(debug_insn 13 12 14 2 sequenceUtil.cpp:185 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 14 13 15 2 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 15 14 16 2 (var_location:DI rev_kmer (const_int 0 [0x0])) -1 (nil))

(debug_insn 16 15 17 2 (var_location:DI kmer (debug_expr:DI D#117)) -1 (nil))

(insn 17 16 18 2 sequenceUtil.cpp:185 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 65 [ kmer_length ])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 18 17 19 2 sequenceUtil.cpp:185 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 23)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9100 [0x238c])
        (nil)))
;; End of basic block 2 -> ( 4 3)

;; Succ edge  4 [91.0%] 
;; Succ edge  3 [9.0%]  (fallthru)

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [9.0%]  (fallthru)
(note 19 18 20 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 20 19 21 3 sequenceUtil.cpp:185 (set (reg/v:DI 62 [ rev_kmer ])
        (const_int 0 [0x0])) -1 (nil))

(jump_insn 21 20 22 3 sequenceUtil.cpp:185 (set (pc)
        (label_ref 46)) -1 (nil))
;; End of basic block 3 -> ( 6)

;; Succ edge  6 [100.0%] 

(barrier 22 21 23)

;; Start of basic block ( 2) -> 4
;; Pred edge  2 [91.0%] 
(code_label 23 22 24 4 19 "" [1 uses])

(note 24 23 25 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 25 24 26 4 sequenceUtil.cpp:184 (set (reg/v:DI 61 [ kmer ])
        (not:DI (reg/v:DI 64 [ kmer_val ]))) -1 (nil))

(insn 26 25 27 4 sequenceUtil.cpp:184 (set (reg/v:SI 60 [ i ])
        (const_int 0 [0x0])) -1 (nil))

(insn 27 26 43 4 sequenceUtil.cpp:184 (set (reg/v:DI 62 [ rev_kmer ])
        (const_int 0 [0x0])) -1 (nil))
;; End of basic block 4 -> ( 5)

;; Succ edge  5 [100.0%]  (fallthru)

;; Start of basic block ( 5 4) -> 5
;; Pred edge  5 [91.0%]  (dfs_back)
;; Pred edge  4 [100.0%]  (fallthru)
(code_label 43 27 28 5 21 "" [1 uses])

(note 28 43 29 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(debug_insn 29 28 30 5 sequenceUtil.cpp:187 (var_location:SI base (and:SI (subreg:SI (reg/v:DI 61 [ kmer ]) 0)
        (const_int 3 [0x3]))) -1 (nil))

(insn 30 29 31 5 sequenceUtil.cpp:188 (parallel [
            (set (reg/v:DI 59 [ rev_kmer ])
                (ashift:DI (reg/v:DI 62 [ rev_kmer ])
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 31 30 32 5 sequenceUtil.cpp:188 (var_location:DI rev_kmer (reg/v:DI 59 [ rev_kmer ])) -1 (nil))

(insn 32 31 33 5 sequenceUtil.cpp:189 (parallel [
            (set (reg:SI 66)
                (and:SI (subreg:SI (reg/v:DI 61 [ kmer ]) 0)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 33 32 34 5 sequenceUtil.cpp:189 (set (reg:DI 67)
        (sign_extend:DI (reg:SI 66))) -1 (nil))

(insn 34 33 35 5 sequenceUtil.cpp:189 (parallel [
            (set (reg/v:DI 62 [ rev_kmer ])
                (plus:DI (reg:DI 67)
                    (reg/v:DI 59 [ rev_kmer ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 35 34 36 5 sequenceUtil.cpp:189 (var_location:DI rev_kmer (reg/v:DI 62 [ rev_kmer ])) -1 (nil))

(insn 36 35 37 5 sequenceUtil.cpp:190 (parallel [
            (set (reg/v:DI 61 [ kmer ])
                (lshiftrt:DI (reg/v:DI 61 [ kmer ])
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 37 36 38 5 sequenceUtil.cpp:190 (var_location:DI kmer (reg/v:DI 61 [ kmer ])) -1 (nil))

(insn 38 37 39 5 sequenceUtil.cpp:185 (parallel [
            (set (reg/v:SI 60 [ i ])
                (plus:SI (reg/v:SI 60 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 39 38 40 5 sequenceUtil.cpp:185 (var_location:SI i (reg/v:SI 60 [ i ])) -1 (nil))

(debug_insn 40 39 41 5 (var_location:SI i (reg/v:SI 60 [ i ])) -1 (nil))

(debug_insn 41 40 42 5 (var_location:DI rev_kmer (reg/v:DI 62 [ rev_kmer ])) -1 (nil))

(debug_insn 42 41 44 5 (var_location:DI kmer (reg/v:DI 61 [ kmer ])) -1 (nil))

(insn 44 42 45 5 sequenceUtil.cpp:185 (set (reg:CC 17 flags)
        (compare:CC (reg/v:SI 65 [ kmer_length ])
            (reg/v:SI 60 [ i ]))) -1 (nil))

(jump_insn 45 44 46 5 sequenceUtil.cpp:185 (set (pc)
        (if_then_else (gtu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 43)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9100 [0x238c])
        (nil)))
;; End of basic block 5 -> ( 5 6)

;; Succ edge  5 [91.0%]  (dfs_back)
;; Succ edge  6 [9.0%]  (fallthru)

;; Start of basic block ( 5 3) -> 6
;; Pred edge  5 [9.0%]  (fallthru)
;; Pred edge  3 [100.0%] 
(code_label 46 45 47 6 20 "" [1 uses])

(note 47 46 48 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(debug_insn 48 47 49 6 sequenceUtil.cpp:398 (var_location:DI rev_kmer (reg/v:DI 62 [ rev_kmer ])) -1 (nil))

(debug_insn 49 48 50 6 (var_location:DI kmer_val (umax:DI (reg/v:DI 62 [ rev_kmer ])
        (reg/v:DI 64 [ kmer_val ]))) -1 (nil))

(insn 50 49 51 6 sequenceUtil.cpp:185 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 62 [ rev_kmer ])
            (reg/v:DI 64 [ kmer_val ]))) -1 (nil))

(insn 51 50 52 6 sequenceUtil.cpp:185 (set (reg:DI 68)
        (if_then_else:DI (geu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (reg/v:DI 62 [ rev_kmer ])
            (reg/v:DI 64 [ kmer_val ]))) -1 (nil))

(insn 52 51 56 6 sequenceUtil.cpp:185 (set (reg:DI 63 [ <result> ])
        (reg:DI 68)) -1 (nil))

(insn 56 52 62 6 sequenceUtil.cpp:405 (set (reg/i:DI 0 ax)
        (reg:DI 63 [ <result> ])) -1 (nil))

(insn 62 56 0 6 sequenceUtil.cpp:405 (use (reg/i:DI 0 ax)) -1 (nil))
;; End of basic block 6 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function (static initializers for sequenceUtil.cpp) (_GLOBAL__I__int_to_base)



try_optimize_cfg iteration 1

Deleting fallthru block 2.
deleting block 2
merging block 4 into block 3
Merged 3 and 4 without moving.
merging block 5 into block 3
Merged 3 and 5 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 2 NOTE_INSN_DELETED)

(note 2 1 4 NOTE_INSN_FUNCTION_BEG)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 4 2 5 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(debug_insn 5 4 6 2 (var_location:SI __initialize_p (const_int 1 [0x1])) -1 (nil))

(debug_insn 6 5 7 2 (var_location:SI __priority (const_int 65535 [0xffff])) -1 (nil))

(insn 7 6 8 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/iostream:72 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZStL8__ioinit") [flags 0x2]  <var_decl 0x7f9012b86e60 __ioinit>)) -1 (nil))

(call_insn 8 7 9 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/iostream:72 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitC1Ev") [flags 0x41]  <function_decl 0x7f901305ea00 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 9 8 10 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/iostream:72 (set (reg:DI 1 dx)
        (symbol_ref:DI ("__dso_handle") [flags 0x40]  <var_decl 0x7f9011dc15a0 __dso_handle>)) -1 (nil))

(insn 10 9 11 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/iostream:72 (set (reg:DI 4 si)
        (symbol_ref:DI ("_ZStL8__ioinit") [flags 0x2]  <var_decl 0x7f9012b86e60 __ioinit>)) -1 (nil))

(insn 11 10 12 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/iostream:72 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZNSt8ios_base4InitD1Ev") [flags 0x41]  <function_decl 0x7f901305ec00 __comp_dtor >)) -1 (nil))

(call_insn 12 11 13 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/iostream:72 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_atexit") [flags 0x41]  <function_decl 0x7f9011dc4200 __cxa_atexit>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(debug_insn 13 12 14 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/iostream:72 (var_location:DI this (symbol_ref:DI ("currAccession") [flags 0x2]  <var_decl 0x7f9012470c80 currAccession>)) -1 (nil))

(debug_insn 14 13 15 2 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 15 14 16 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/allocator.h:101 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 16 15 17 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:178 (var_location:DI __p (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>)) -1 (nil))

(debug_insn 17 16 18 2 (var_location:DI this (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>)) -1 (nil))

(debug_insn 18 17 19 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2147 (var_location:DI D.4294967278 (symbol_ref:DI ("currAccession") [flags 0x2]  <var_decl 0x7f9012470c80 currAccession>)) -1 (nil))

(debug_insn 19 18 20 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2147 (var_location:DI this (debug_expr:DI D#18)) -1 (nil))

(debug_insn 20 19 21 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2147 (var_location:DI __dat (plus:DI (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>)
        (const_int 24 [0x18]))) -1 (nil))

(debug_insn 21 20 22 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2147 (var_location:DI __a (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 22 21 23 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:260 (var_location:DI D.4294967282 (debug_expr:DI D#18)) -1 (nil))

(debug_insn 23 22 24 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:260 (var_location:DI this (debug_expr:DI D#14)) -1 (nil))

(debug_insn 24 23 25 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:260 (var_location:DI __a (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 25 24 26 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/allocator.h:104 (var_location:DI this (debug_expr:DI D#14)) -1 (nil))

(debug_insn 26 25 27 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/allocator.h:104 (var_location:DI D.39176 (clobber (const_int 0 [0x0]))) -1 (nil))

(insn 27 26 28 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:260 (set (reg/f:DI 58)
        (symbol_ref:DI ("currAccession") [flags 0x2]  <var_decl 0x7f9012470c80 currAccession>)) -1 (nil))

(insn 28 27 32 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:260 (set (mem/s/f/c:DI (reg/f:DI 58) [32 currAccession._M_dataplus._M_p+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>)
                (const_int 24 [0x18])))) -1 (nil))

(debug_insn 32 28 33 2 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 33 32 34 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/allocator.h:109 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(insn 34 33 35 2 sequenceUtil.cpp:28 (set (reg:DI 1 dx)
        (symbol_ref:DI ("__dso_handle") [flags 0x40]  <var_decl 0x7f9011dc15a0 __dso_handle>)) -1 (nil))

(insn 35 34 36 2 sequenceUtil.cpp:28 (set (reg:DI 4 si)
        (symbol_ref:DI ("currAccession") [flags 0x2]  <var_decl 0x7f9012470c80 currAccession>)) -1 (nil))

(insn 36 35 37 2 sequenceUtil.cpp:28 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZNSsD1Ev") [flags 0x41]  <function_decl 0x7f901329ba00 __comp_dtor >)) -1 (nil))

(call_insn/j 37 36 38 2 sequenceUtil.cpp:28 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_atexit") [flags 0x41]  <function_decl 0x7f9011dc4200 __cxa_atexit>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 2 -> ( 1)

;; Succ edge  EXIT [100.0%]  (ab,sibcall)

(barrier 38 37 0)


;; Function static _CharT* std::basic_string<_CharT, _Traits, _Alloc>::_S_construct(_InIterator, _InIterator, const _Alloc&, std::forward_iterator_tag) [with _FwdIterator = char*, _CharT = char, _Traits = std::char_traits<char>, _Alloc = std::allocator<char>] (_ZNSs12_S_constructIPcEES0_T_S1_RKSaIcESt20forward_iterator_tag)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
Forwarding edge 2->4 to 5 failed.
merging block 4 into block 2
Merged 2 and 4 without moving.
merging block 5 into block 2
Merged 2 and 5 without moving.
Forwarding edge 15->16 to 11 failed.
merging block 19 into block 18
Merged 18 and 19 without moving.


try_optimize_cfg iteration 2

Forwarding edge 15->16 to 11 failed.


try_optimize_cfg iteration 1

Forwarding edge 12->13 to 8 failed.
(note 1 0 7 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 7 1 2 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 2 7 3 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:123 (set (reg/v/f:DI 62 [ __beg ])
        (reg:DI 5 di [ __beg ])) -1 (nil))

(insn 3 2 4 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:123 (set (reg/v/f:DI 63 [ __end ])
        (reg:DI 4 si [ __end ])) -1 (nil))

(insn 4 3 5 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:123 (set (reg/v/f:DI 64 [ __a ])
        (reg:DI 1 dx [ __a ])) -1 (nil))

(insn 5 4 6 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:123 (set (reg/v:QI 65 [ D.37805 ])
        (mem/s/c:QI (reg/f:DI 53 virtual-incoming-args) [86 D.37805+0 S1 A64])) -1 (expr_list:REG_EQUIV (mem/s/c:QI (reg/f:DI 53 virtual-incoming-args) [86 D.37805+0 S1 A64])
        (nil)))

(note 6 5 12 2 NOTE_INSN_FUNCTION_BEG)

(insn 12 6 13 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:128 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 62 [ __beg ])
            (reg/v/f:DI 63 [ __end ]))) -1 (nil))

(jump_insn 13 12 14 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:128 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 20)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 6100 [0x17d4])
        (nil)))
;; End of basic block 2 -> ( 3 4)

;; Succ edge  3 [39.0%]  (fallthru)
;; Succ edge  4 [61.0%] 

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [39.0%]  (fallthru)
(note 14 13 15 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(debug_insn 15 14 16 3 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:178 (var_location:DI __p (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>)) -1 (nil))

(debug_insn 16 15 17 3 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:129 (var_location:DI this (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>)) -1 (nil))

(insn 17 16 18 3 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (set (reg/f:DI 58 [ D.39749 ])
        (const:DI (plus:DI (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>)
                (const_int 24 [0x18])))) -1 (nil))

(jump_insn 18 17 19 3 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:129 (set (pc)
        (label_ref 67)) -1 (nil))
;; End of basic block 3 -> ( 11)

;; Succ edge  11 [100.0%] 

(barrier 19 18 20)

;; Start of basic block ( 2) -> 4
;; Pred edge  2 [61.0%] 
(code_label 20 19 21 4 31 "" [1 uses])

(note 21 20 22 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(debug_insn 22 21 23 4 (var_location:DI __ptr (reg/v/f:DI 62 [ __beg ])) -1 (nil))

(insn 23 22 24 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 62 [ __beg ])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 24 23 25 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:132 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 72)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 8500 [0x2134])
        (nil)))
;; End of basic block 4 -> ( 5 12)

;; Succ edge  5 [15.0%]  (fallthru)
;; Succ edge  12 [85.0%] 

;; Start of basic block ( 4) -> 5
;; Pred edge  4 [15.0%]  (fallthru)
(note 25 24 26 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 26 25 27 5 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:132 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 63 [ __end ])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 27 26 28 5 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:132 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 102)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 1500 [0x5dc])
        (nil)))
;; End of basic block 5 -> ( 6 14)

;; Succ edge  6 [85.0%]  (fallthru)
;; Succ edge  14 [15.0%] 

;; Start of basic block ( 5) -> 6
;; Pred edge  5 [85.0%]  (fallthru)
(note 28 27 29 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 29 28 30 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:134 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC0") [flags 0x2]  <string_cst 0x7f901197d730>)) -1 (nil))

(call_insn 30 29 31 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:134 (call (mem:QI (symbol_ref:DI ("_ZSt19__throw_logic_errorPKc") [flags 0x41]  <function_decl 0x7f9013bc4500 __throw_logic_error>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 6 -> ()


(barrier 31 30 95)

;; Start of basic block ( 12) -> 7
;; Pred edge  12 [38.8%] 
(code_label 95 31 32 7 36 "" [1 uses])

(note 32 95 33 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(debug_insn 33 32 34 7 (var_location:DI __c1 (reg/f:DI 58 [ D.39749 ])) -1 (nil))

(debug_insn 34 33 35 7 (var_location:DI __c2 (reg/v/f:DI 62 [ __beg ])) -1 (nil))

(insn 35 34 36 7 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:247 (set (reg:QI 66)
        (mem:QI (reg/v/f:DI 62 [ __beg ]) [0 S1 A8])) -1 (nil))

(insn 36 35 37 7 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:247 (set (mem:QI (reg/f:DI 58 [ D.39749 ]) [0 S1 A8])
        (reg:QI 66)) -1 (nil))

(jump_insn 37 36 38 7 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:247 (set (pc)
        (label_ref 52)) -1 (nil))
;; End of basic block 7 -> ( 9)

;; Succ edge  9 [100.0%] 

(barrier 38 37 98)

;; Start of basic block ( 14 13) -> 8
;; Pred edge  14 [100.0%] 
;; Pred edge  13 [100.0%] 
(code_label 98 38 39 8 37 "" [2 uses])

(note 39 98 40 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(debug_insn 40 39 41 8 (var_location:DI __s1 (reg/f:DI 58 [ D.39749 ])) -1 (nil))

(debug_insn 41 40 42 8 (var_location:DI __s2 (reg/v/f:DI 62 [ __beg ])) -1 (nil))

(debug_insn 42 41 43 8 (var_location:DI __n (reg/v:DI 60 [ __dnew ])) -1 (nil))

(insn 43 42 44 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 67)
        (reg/f:DI 58 [ D.39749 ])) -1 (nil))

(insn 44 43 45 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 68)
        (reg/v/f:DI 62 [ __beg ])) -1 (nil))

(insn 45 44 46 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 69)
        (reg/v:DI 60 [ __dnew ])) -1 (nil))

(insn 46 45 47 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 1 dx)
        (reg:DI 69)) -1 (nil))

(insn 47 46 48 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 4 si)
        (reg:DI 68)) -1 (nil))

(insn 48 47 49 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 5 di)
        (reg:DI 67)) -1 (nil))

(call_insn 49 48 50 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memcpy") [flags 0x41]  <function_decl 0x7f9011be1f00 memcpy>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 50 49 51 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 70)
        (reg:DI 0 ax)) -1 (nil))

(insn 51 50 52 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:275 (set (reg:DI 71)
        (reg:DI 70)) -1 (nil))
;; End of basic block 8 -> ( 9)

;; Succ edge  9 [100.0%]  (fallthru)

;; Start of basic block ( 8 7) -> 9
;; Pred edge  8 [100.0%]  (fallthru)
;; Pred edge  7 [100.0%] 
(code_label 52 51 53 9 35 "" [1 uses])

(note 53 52 54 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(debug_insn 54 53 55 9 (var_location:DI this (reg/v/f:DI 59 [ __r ])) -1 (nil))

(debug_insn 55 54 56 9 (var_location:DI __n (reg/v:DI 60 [ __dnew ])) -1 (nil))

(debug_insn 56 55 57 9 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:178 (var_location:DI __p (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>)) -1 (nil))

(insn 57 56 58 9 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:202 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 59 [ __r ])
            (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40]  <var_decl 0x7f90132c3280 _S_empty_rep_storage>))) -1 (nil))

(jump_insn 58 57 59 9 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:202 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 67)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9996 [0x270c])
        (nil)))
;; End of basic block 9 -> ( 10 11)

;; Succ edge  10 [0.0%]  (fallthru)
;; Succ edge  11 [100.0%] 

;; Start of basic block ( 9) -> 10
;; Pred edge  9 [0.0%]  (fallthru)
(note 59 58 60 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(debug_insn 60 59 61 10 (var_location:DI this (reg/v/f:DI 59 [ __r ])) -1 (nil))

(insn 61 60 62 10 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:196 (set (mem/s:SI (plus:DI (reg/v/f:DI 59 [ __r ])
                (const_int 16 [0x10])) [5 <variable>.D.11486._M_refcount+0 S4 A64])
        (const_int 0 [0x0])) -1 (nil))

(insn 62 61 63 10 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:206 (set (mem/s:DI (reg/v/f:DI 59 [ __r ]) [14 <variable>.D.11486._M_length+0 S8 A64])
        (reg/v:DI 60 [ __dnew ])) -1 (nil))

(debug_insn 63 62 64 10 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:206 (var_location:DI this (reg/v/f:DI 59 [ __r ])) -1 (nil))

(debug_insn 64 63 65 10 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:207 (var_location:DI __c1 (plus:DI (reg/f:DI 58 [ D.39749 ])
        (reg/v:DI 60 [ __dnew ]))) -1 (nil))

(debug_insn 65 64 66 10 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:207 (var_location:DI __c2 (clobber (const_int 0 [0x0]))) -1 (nil))

(insn 66 65 67 10 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:247 (set (mem:QI (plus:DI (reg/f:DI 58 [ D.39749 ])
                (reg/v:DI 60 [ __dnew ])) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil))
;; End of basic block 10 -> ( 11)

;; Succ edge  11 [100.0%]  (fallthru)

;; Start of basic block ( 3 10 9) -> 11
;; Pred edge  3 [100.0%] 
;; Pred edge  10 [100.0%]  (fallthru)
;; Pred edge  9 [100.0%] 
(code_label 67 66 68 11 32 "" [2 uses])

(note 68 67 69 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 69 68 70 11 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:247 (set (reg:DI 61 [ <result> ])
        (reg/f:DI 58 [ D.39749 ])) -1 (nil))

(jump_insn 70 69 71 11 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:247 (set (pc)
        (label_ref 127)) -1 (nil))
;; End of basic block 11 -> ( 15)

;; Succ edge  15 [100.0%] 

(barrier 71 70 72)

;; Start of basic block ( 4) -> 12
;; Pred edge  4 [85.0%] 
(code_label 72 71 73 12 33 "" [1 uses])

(note 73 72 74 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(debug_insn 74 73 75 12 (var_location:DI __first (reg/v/f:DI 62 [ __beg ])) -1 (nil))

(debug_insn 75 74 76 12 (var_location:DI __last (reg/v/f:DI 63 [ __end ])) -1 (nil))

(debug_insn 76 75 77 12 (var_location:DI D.39759 (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 77 76 78 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_iterator_base_funcs.h:114 (var_location:DI __first (reg/v/f:DI 62 [ __beg ])) -1 (nil))

(debug_insn 78 77 79 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_iterator_base_funcs.h:114 (var_location:DI __last (reg/v/f:DI 63 [ __end ])) -1 (nil))

(insn 79 78 80 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:137 (parallel [
            (set (reg/v:DI 60 [ __dnew ])
                (minus:DI (reg/v/f:DI 63 [ __end ])
                    (reg/v/f:DI 62 [ __beg ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 80 79 81 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:137 (var_location:DI __dnew (reg/v:DI 60 [ __dnew ])) -1 (nil))

(insn 81 80 82 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg:DI 1 dx)
        (reg/v/f:DI 64 [ __a ])) -1 (nil))

(insn 82 81 83 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg:DI 4 si)
        (const_int 0 [0x0])) -1 (nil))

(insn 83 82 84 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg:DI 5 di)
        (reg/v:DI 60 [ __dnew ])) -1 (nil))

(call_insn 84 83 85 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep9_S_createEmmRKSaIcE") [flags 0x41]  <function_decl 0x7f90132c9400 _S_create>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 85 84 86 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg/v/f:DI 59 [ __r ])
        (reg:DI 0 ax)) -1 (nil))

(debug_insn 86 85 87 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (var_location:DI __r (reg/v/f:DI 59 [ __r ])) -1 (nil))

(debug_insn 87 86 88 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (var_location:DI this (reg/v/f:DI 59 [ __r ])) -1 (nil))

(insn 88 87 89 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (parallel [
            (set (reg/f:DI 58 [ D.39749 ])
                (plus:DI (reg/v/f:DI 59 [ __r ])
                    (const_int 24 [0x18])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 89 88 90 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (var_location:DI __p (reg/f:DI 58 [ D.39749 ])) -1 (nil))

(debug_insn 90 89 91 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (var_location:DI __k1 (reg/v/f:DI 62 [ __beg ])) -1 (nil))

(debug_insn 91 90 92 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (var_location:DI __k2 (reg/v/f:DI 63 [ __end ])) -1 (nil))

(debug_insn 92 91 93 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:385 (var_location:DI __d (reg/f:DI 58 [ D.39749 ])) -1 (nil))

(debug_insn 93 92 94 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:385 (var_location:DI __s (reg/v/f:DI 62 [ __beg ])) -1 (nil))

(debug_insn 94 93 96 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:385 (var_location:DI __n (reg/v:DI 60 [ __dnew ])) -1 (nil))

(insn 96 94 97 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:341 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:DI 60 [ __dnew ])
            (const_int 1 [0x1]))) -1 (nil))

(jump_insn 97 96 101 12 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:341 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 95)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 3884 [0xf2c])
        (nil)))
;; End of basic block 12 -> ( 7 13)

;; Succ edge  7 [38.8%] 
;; Succ edge  13 [61.2%]  (fallthru)

;; Start of basic block ( 12) -> 13
;; Pred edge  12 [61.2%]  (fallthru)
(note 101 97 99 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(jump_insn 99 101 100 13 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:341 (set (pc)
        (label_ref 98)) -1 (nil))
;; End of basic block 13 -> ( 8)

;; Succ edge  8 [100.0%] 

(barrier 100 99 102)

;; Start of basic block ( 5) -> 14
;; Pred edge  5 [15.0%] 
(code_label 102 100 103 14 34 "" [1 uses])

(note 103 102 104 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(debug_insn 104 103 105 14 (var_location:DI __first (const_int 0 [0x0])) -1 (nil))

(debug_insn 105 104 106 14 (var_location:DI __last (const_int 0 [0x0])) -1 (nil))

(debug_insn 106 105 107 14 (var_location:DI D.39759 (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 107 106 108 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_iterator_base_funcs.h:114 (var_location:DI __first (const_int 0 [0x0])) -1 (nil))

(debug_insn 108 107 109 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_iterator_base_funcs.h:114 (var_location:DI __last (const_int 0 [0x0])) -1 (nil))

(debug_insn 109 108 110 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:137 (var_location:DI __dnew (const_int 0 [0x0])) -1 (nil))

(insn 110 109 111 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg:DI 1 dx)
        (reg/v/f:DI 64 [ __a ])) -1 (nil))

(insn 111 110 112 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg:DI 4 si)
        (const_int 0 [0x0])) -1 (nil))

(insn 112 111 113 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg:DI 5 di)
        (const_int 0 [0x0])) -1 (nil))

(call_insn 113 112 114 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep9_S_createEmmRKSaIcE") [flags 0x41]  <function_decl 0x7f90132c9400 _S_create>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 114 113 115 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (set (reg/v/f:DI 59 [ __r ])
        (reg:DI 0 ax)) -1 (nil))

(debug_insn 115 114 116 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (var_location:DI __r (reg/v/f:DI 59 [ __r ])) -1 (nil))

(debug_insn 116 115 117 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:139 (var_location:DI this (reg/v/f:DI 59 [ __r ])) -1 (nil))

(insn 117 116 118 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (parallel [
            (set (reg/f:DI 58 [ D.39749 ])
                (plus:DI (reg/v/f:DI 59 [ __r ])
                    (const_int 24 [0x18])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 118 117 119 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (var_location:DI __p (reg/f:DI 58 [ D.39749 ])) -1 (nil))

(debug_insn 119 118 120 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (var_location:DI __k1 (const_int 0 [0x0])) -1 (nil))

(debug_insn 120 119 121 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (var_location:DI __k2 (const_int 0 [0x0])) -1 (nil))

(debug_insn 121 120 122 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:385 (var_location:DI __d (reg/f:DI 58 [ D.39749 ])) -1 (nil))

(debug_insn 122 121 123 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:385 (var_location:DI __s (const_int 0 [0x0])) -1 (nil))

(debug_insn 123 122 124 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:385 (var_location:DI __n (const_int 0 [0x0])) -1 (nil))

(insn 124 123 125 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (set (reg/v:DI 60 [ __dnew ])
        (const_int 0 [0x0])) -1 (nil))

(jump_insn 125 124 126 14 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:215 (set (pc)
        (label_ref 98)) -1 (nil))
;; End of basic block 14 -> ( 8)

;; Succ edge  8 [100.0%] 

(barrier 126 125 127)

;; Start of basic block ( 11) -> 15
;; Pred edge  11 [100.0%] 
(code_label 127 126 137 15 28 "" [1 uses])

(note 137 127 128 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(insn 128 137 134 15 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:149 (set (reg/i:DI 0 ax)
        (reg:DI 61 [ <result> ])) -1 (nil))

(insn 134 128 0 15 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.tcc:149 (use (reg/i:DI 0 ax)) -1 (nil))
;; End of basic block 15 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function float compute_entropy(kmer_int_type_t, unsigned int) (_Z15compute_entropyyj)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
Edge 14->16 redirected to 17
Forwarding edge 14->15 to 18 failed.
Forwarding edge 14->15 to 18 failed.
deleting block 16
Edge 18->20 redirected to 21
Forwarding edge 18->19 to 22 failed.
Forwarding edge 18->19 to 22 failed.
deleting block 20
Edge 22->24 redirected to 25
Forwarding edge 22->23 to 26 failed.
Forwarding edge 22->23 to 26 failed.
deleting block 24
merging block 27 into block 26
Merged 26 and 27 without moving.
merging block 28 into block 26
Merged 26 and 28 without moving.


try_optimize_cfg iteration 2

Forwarding edge 14->15 to 18 failed.
Forwarding edge 18->19 to 22 failed.
Forwarding edge 22->23 to 26 failed.


try_optimize_cfg iteration 1

Forwarding edge 13->14 to 16 failed.
Forwarding edge 16->17 to 19 failed.
Forwarding edge 19->20 to 22 failed.
(note 2 0 6 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 6 2 3 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 2 sequenceUtil.cpp:336 (set (reg/v:DI 73 [ kmer ])
        (reg:DI 5 di [ kmer ])) -1 (nil))

(insn 4 3 5 2 sequenceUtil.cpp:336 (set (reg/v:SI 74 [ kmer_length ])
        (reg:SI 4 si [ kmer_length ])) -1 (nil))

(note 5 4 8 2 NOTE_INSN_FUNCTION_BEG)

(insn 8 5 9 2 sequenceUtil.cpp:338 (set (mem/s/j:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -16 [0xfffffffffffffff0])) [0 counts+0 S1 A128])
        (const_int 0 [0x0])) -1 (nil))

(insn 9 8 10 2 sequenceUtil.cpp:338 (set (mem/s/j:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -15 [0xfffffffffffffff1])) [0 counts+1 S1 A8])
        (const_int 0 [0x0])) -1 (nil))

(insn 10 9 11 2 sequenceUtil.cpp:338 (set (mem/s/j:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -14 [0xfffffffffffffff2])) [0 counts+2 S1 A16])
        (const_int 0 [0x0])) -1 (nil))

(insn 11 10 12 2 sequenceUtil.cpp:338 (set (mem/s/j:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -13 [0xfffffffffffffff3])) [0 counts+3 S1 A8])
        (const_int 0 [0x0])) -1 (nil))

(debug_insn 12 11 13 2 sequenceUtil.cpp:340 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 13 12 14 2 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 14 13 15 2 (var_location:DI kmer (reg/v:DI 73 [ kmer ])) -1 (nil))

(insn 15 14 16 2 sequenceUtil.cpp:340 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 74 [ kmer_length ])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 16 15 17 2 sequenceUtil.cpp:340 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 21)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9100 [0x238c])
        (nil)))
;; End of basic block 2 -> ( 4 3)

;; Succ edge  4 [91.0%] 
;; Succ edge  3 [9.0%]  (fallthru)

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [9.0%]  (fallthru)
(note 17 16 18 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 18 17 19 3 sequenceUtil.cpp:340 (set (reg:QI 65 [ prephitmp.698 ])
        (const_int 0 [0x0])) -1 (nil))

(jump_insn 19 18 20 3 sequenceUtil.cpp:340 (set (pc)
        (label_ref 43)) -1 (nil))
;; End of basic block 3 -> ( 7)

;; Succ edge  7 [100.0%] 

(barrier 20 19 21)

;; Start of basic block ( 2) -> 4
;; Pred edge  2 [91.0%] 
(code_label 21 20 22 4 43 "" [1 uses])

(note 22 21 23 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 23 22 38 4 sequenceUtil.cpp:340 (set (reg/v:SI 68 [ i ])
        (const_int 0 [0x0])) -1 (nil))
;; End of basic block 4 -> ( 5)

;; Succ edge  5 [100.0%]  (fallthru)

;; Start of basic block ( 4 5) -> 5
;; Pred edge  4 [100.0%]  (fallthru)
;; Pred edge  5 [91.0%]  (dfs_back)
(code_label 38 23 24 5 45 "" [1 uses])

(note 24 38 25 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 25 24 26 5 sequenceUtil.cpp:342 (parallel [
            (set (reg/v:SI 67 [ c ])
                (and:SI (subreg:SI (reg/v:DI 73 [ kmer ]) 0)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 26 25 27 5 sequenceUtil.cpp:342 (var_location:SI c (reg/v:SI 67 [ c ])) -1 (nil))

(insn 27 26 28 5 sequenceUtil.cpp:343 (parallel [
            (set (reg/v:DI 73 [ kmer ])
                (lshiftrt:DI (reg/v:DI 73 [ kmer ])
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 28 27 29 5 sequenceUtil.cpp:343 (var_location:DI kmer (reg/v:DI 73 [ kmer ])) -1 (nil))

(insn 29 28 30 5 sequenceUtil.cpp:344 (set (reg:DI 75)
        (sign_extend:DI (reg/v:SI 67 [ c ]))) -1 (nil))

(insn 30 29 31 5 sequenceUtil.cpp:344 (set (reg:DI 76)
        (sign_extend:DI (reg/v:SI 67 [ c ]))) -1 (nil))

(insn 31 30 32 5 sequenceUtil.cpp:344 (set (reg:QI 78)
        (mem/s/j:QI (plus:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (reg:DI 76))
                (const_int -16 [0xfffffffffffffff0])) [0 counts S1 A8])) -1 (nil))

(insn 32 31 33 5 sequenceUtil.cpp:344 (parallel [
            (set (reg:QI 77)
                (plus:QI (reg:QI 78)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 33 32 34 5 sequenceUtil.cpp:344 (set (mem/s/j:QI (plus:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (reg:DI 75))
                (const_int -16 [0xfffffffffffffff0])) [0 counts S1 A8])
        (reg:QI 77)) -1 (expr_list:REG_EQUAL (plus:QI (mem/s/j:QI (plus:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                        (reg:DI 76))
                    (const_int -16 [0xfffffffffffffff0])) [0 counts S1 A8])
            (const_int 1 [0x1]))
        (nil)))

(insn 34 33 35 5 sequenceUtil.cpp:340 (parallel [
            (set (reg/v:SI 68 [ i ])
                (plus:SI (reg/v:SI 68 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 35 34 36 5 sequenceUtil.cpp:340 (var_location:SI i (reg/v:SI 68 [ i ])) -1 (nil))

(debug_insn 36 35 37 5 (var_location:SI i (reg/v:SI 68 [ i ])) -1 (nil))

(debug_insn 37 36 39 5 (var_location:DI kmer (reg/v:DI 73 [ kmer ])) -1 (nil))

(insn 39 37 40 5 sequenceUtil.cpp:340 (set (reg:CC 17 flags)
        (compare:CC (reg/v:SI 74 [ kmer_length ])
            (reg/v:SI 68 [ i ]))) -1 (nil))

(jump_insn 40 39 41 5 sequenceUtil.cpp:340 (set (pc)
        (if_then_else (gtu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 38)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9100 [0x238c])
        (nil)))
;; End of basic block 5 -> ( 5 6)

;; Succ edge  5 [91.0%]  (dfs_back)
;; Succ edge  6 [9.0%]  (fallthru)

;; Start of basic block ( 5) -> 6
;; Pred edge  5 [9.0%]  (fallthru)
(note 41 40 42 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 42 41 43 6 sequenceUtil.cpp:340 (set (reg:QI 65 [ prephitmp.698 ])
        (mem/s/j:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -16 [0xfffffffffffffff0])) [0 counts+0 S1 A128])) -1 (nil))
;; End of basic block 6 -> ( 7)

;; Succ edge  7 [100.0%]  (fallthru)

;; Start of basic block ( 6 3) -> 7
;; Pred edge  6 [100.0%]  (fallthru)
;; Pred edge  3 [100.0%] 
(code_label 43 42 44 7 44 "" [1 uses])

(note 44 43 45 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(debug_insn 45 44 46 7 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 46 45 47 7 (var_location:SF entropy (const_double:SF 0.0 [0x0.0p+0])) -1 (nil))

(insn 47 46 48 7 sequenceUtil.cpp:351 (set (reg:DI 79)
        (zero_extend:DI (reg/v:SI 74 [ kmer_length ]))) -1 (nil))

(insn 48 47 49 7 sequenceUtil.cpp:351 (set (reg:CCGOC 17 flags)
        (compare:CCGOC (reg:DI 79)
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 49 48 199 7 sequenceUtil.cpp:351 (set (pc)
        (if_then_else (lt (reg:CCGOC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 53)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 2100 [0x834])
        (nil)))
;; End of basic block 7 -> ( 9 8)

;; Succ edge  9 [21.0%] 
;; Succ edge  8 [79.0%]  (fallthru)

;; Start of basic block ( 7) -> 8
;; Pred edge  7 [79.0%]  (fallthru)
(note 199 49 50 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 50 199 51 8 sequenceUtil.cpp:351 (set (reg:SF 71 [ D.32917 ])
        (float:SF (reg:DI 79))) -1 (nil))

(jump_insn 51 50 52 8 sequenceUtil.cpp:351 (set (pc)
        (label_ref 59)) -1 (nil))
;; End of basic block 8 -> ( 10)

;; Succ edge  10 [100.0%] 

(barrier 52 51 53)

;; Start of basic block ( 7) -> 9
;; Pred edge  7 [21.0%] 
(code_label 53 52 200 9 46 "" [1 uses])

(note 200 53 54 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 54 200 55 9 sequenceUtil.cpp:351 (parallel [
            (set (reg:DI 81)
                (lshiftrt:DI (reg:DI 79)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 55 54 56 9 sequenceUtil.cpp:351 (parallel [
            (set (reg:DI 82)
                (and:DI (reg:DI 79)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 56 55 57 9 sequenceUtil.cpp:351 (parallel [
            (set (reg:DI 81)
                (ior:DI (reg:DI 81)
                    (reg:DI 82)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 57 56 58 9 sequenceUtil.cpp:351 (set (reg:SF 80)
        (float:SF (reg:DI 81))) -1 (nil))

(insn 58 57 59 9 sequenceUtil.cpp:351 (set (reg:SF 71 [ D.32917 ])
        (plus:SF (reg:SF 80)
            (reg:SF 80))) -1 (nil))
;; End of basic block 9 -> ( 10)

;; Succ edge  10 [100.0%]  (fallthru)

;; Start of basic block ( 8 9) -> 10
;; Pred edge  8 [100.0%] 
;; Pred edge  9 [100.0%]  (fallthru)
(code_label 59 58 201 10 47 "" [1 uses])

(note 201 59 60 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(insn 60 201 61 10 sequenceUtil.cpp:351 (set (reg:SI 84)
        (sign_extend:SI (reg:QI 65 [ prephitmp.698 ]))) -1 (nil))

(insn 61 60 62 10 sequenceUtil.cpp:351 (set (reg:SF 83)
        (float:SF (reg:SI 84))) -1 (nil))

(insn 62 61 63 10 sequenceUtil.cpp:351 (set (reg/v:SF 64 [ prob.718 ])
        (div:SF (reg:SF 83)
            (reg:SF 71 [ D.32917 ]))) -1 (nil))

(debug_insn 63 62 64 10 sequenceUtil.cpp:351 (var_location:SF prob (reg/v:SF 64 [ prob.718 ])) -1 (nil))

(insn 64 63 65 10 sequenceUtil.cpp:353 (set (reg:SF 85)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 0.0 [0x0.0p+0])
        (nil)))

(insn 65 64 66 10 sequenceUtil.cpp:353 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:SF 64 [ prob.718 ])
            (reg:SF 85))) -1 (nil))

(jump_insn 66 65 67 10 sequenceUtil.cpp:353 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 2900 [0xb54])
        (nil)))
;; End of basic block 10 -> ( 12 11)

;; Succ edge  12 [29.0%] 
;; Succ edge  11 [71.0%]  (fallthru)

;; Start of basic block ( 10) -> 11
;; Pred edge  10 [71.0%]  (fallthru)
(note 67 66 68 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 68 67 69 11 sequenceUtil.cpp:353 (set (reg/v:SF 69 [ entropy ])
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 0.0 [0x0.0p+0])
        (nil)))

(jump_insn 69 68 70 11 sequenceUtil.cpp:353 (set (pc)
        (label_ref 88)) -1 (nil))
;; End of basic block 11 -> ( 13)

;; Succ edge  13 [100.0%] 

(barrier 70 69 71)

;; Start of basic block ( 10) -> 12
;; Pred edge  10 [29.0%] 
(code_label 71 70 72 12 48 "" [1 uses])

(note 72 71 73 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(insn 73 72 74 12 sequenceUtil.cpp:354 (set (reg:SF 87)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC2") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 74 73 75 12 sequenceUtil.cpp:354 (set (reg:SF 86)
        (div:SF (reg:SF 87)
            (reg/v:SF 64 [ prob.718 ]))) -1 (nil))

(insn 75 74 76 12 sequenceUtil.cpp:354 (set (reg:DF 88)
        (float_extend:DF (reg:SF 86))) -1 (nil))

(insn 76 75 77 12 sequenceUtil.cpp:354 (set (reg:DF 21 xmm0)
        (reg:DF 88)) -1 (nil))

(call_insn 77 76 78 12 sequenceUtil.cpp:354 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("log") [flags 0x41]  <function_decl 0x7f9013e35c00 log>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
        (nil)))

(insn 78 77 79 12 sequenceUtil.cpp:354 (set (reg:DF 63 [ temp.723 ])
        (reg:DF 21 xmm0)) -1 (nil))

(debug_insn 79 78 80 12 sequenceUtil.cpp:354 (var_location:SF val (float_truncate:SF (div:DF (mult:DF (float_extend:DF (reg/v:SF 64 [ prob.718 ]))
                (reg:DF 63 [ temp.723 ]))
            (const_double:DF 6.9314718055994528622676398299518041312694549560546875e-1 [0x0.b17217f7d1cf78p+0])))) -1 (nil))

(insn 80 79 81 12 sequenceUtil.cpp:355 (set (reg:DF 89)
        (float_extend:DF (reg/v:SF 64 [ prob.718 ]))) -1 (nil))

(insn 81 80 82 12 sequenceUtil.cpp:355 (set (reg:DF 90)
        (mult:DF (reg:DF 89)
            (reg:DF 63 [ temp.723 ]))) -1 (nil))

(insn 82 81 83 12 sequenceUtil.cpp:355 (set (reg:DF 92)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC3") [flags 0x2]) [96 S8 A64])) -1 (expr_list:REG_EQUAL (const_double:DF 6.9314718055994528622676398299518041312694549560546875e-1 [0x0.b17217f7d1cf78p+0])
        (nil)))

(insn 83 82 84 12 sequenceUtil.cpp:355 (set (reg:DF 91)
        (div:DF (reg:DF 90)
            (reg:DF 92))) -1 (nil))

(insn 84 83 85 12 sequenceUtil.cpp:355 (set (reg:SF 93)
        (float_truncate:SF (reg:DF 91))) -1 (nil))

(insn 85 84 86 12 sequenceUtil.cpp:355 (set (reg:SF 94)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 0.0 [0x0.0p+0])
        (nil)))

(insn 86 85 87 12 sequenceUtil.cpp:355 (set (reg/v:SF 69 [ entropy ])
        (plus:SF (reg:SF 93)
            (reg:SF 94))) -1 (nil))

(debug_insn 87 86 88 12 sequenceUtil.cpp:355 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))
;; End of basic block 12 -> ( 13)

;; Succ edge  13 [100.0%]  (fallthru)

;; Start of basic block ( 11 12) -> 13
;; Pred edge  11 [100.0%] 
;; Pred edge  12 [100.0%]  (fallthru)
(code_label 88 87 89 13 49 "" [1 uses])

(note 89 88 90 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(debug_insn 90 89 91 13 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))

(debug_insn 91 90 92 13 sequenceUtil.cpp:349 (var_location:SI i (const_int 1 [0x1])) -1 (nil))

(debug_insn 92 91 93 13 (var_location:SI i (const_int 1 [0x1])) -1 (nil))

(debug_insn 93 92 94 13 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))

(insn 94 93 95 13 sequenceUtil.cpp:351 (set (reg:SI 96)
        (sign_extend:SI (mem/s/j:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -15 [0xfffffffffffffff1])) [0 counts+1 S1 A8]))) -1 (nil))

(insn 95 94 96 13 sequenceUtil.cpp:351 (set (reg:SF 95)
        (float:SF (reg:SI 96))) -1 (nil))

(insn 96 95 97 13 sequenceUtil.cpp:351 (set (reg/v:SF 62 [ prob.728 ])
        (div:SF (reg:SF 95)
            (reg:SF 71 [ D.32917 ]))) -1 (nil))

(debug_insn 97 96 98 13 sequenceUtil.cpp:351 (var_location:SF prob (reg/v:SF 62 [ prob.728 ])) -1 (nil))

(insn 98 97 99 13 sequenceUtil.cpp:353 (set (reg:SF 97)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 0.0 [0x0.0p+0])
        (nil)))

(insn 99 98 100 13 sequenceUtil.cpp:353 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:SF 62 [ prob.728 ])
            (reg:SF 97))) -1 (nil))

(jump_insn 100 99 202 13 sequenceUtil.cpp:353 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 210)
            (pc))) 617 {*jcc_1} (expr_list:REG_BR_PROB (const_int 2900 [0xb54])
        (nil)))
;; End of basic block 13 -> ( 15 14)

;; Succ edge  15 [29.0%] 
;; Succ edge  14 [71.0%]  (fallthru)

;; Start of basic block ( 13) -> 14
;; Pred edge  13 [71.0%]  (fallthru)
(note 202 100 101 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(jump_insn 101 202 102 14 sequenceUtil.cpp:353 (set (pc)
        (label_ref 119)) -1 (nil))
;; End of basic block 14 -> ( 16)

;; Succ edge  16 [100.0%] 

(barrier 102 101 210)

;; Start of basic block ( 13) -> 15
;; Pred edge  13 [29.0%] 
(code_label 210 102 104 15 57 "" [1 uses])

(note 104 210 105 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(insn 105 104 106 15 sequenceUtil.cpp:354 (set (reg:SF 99)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC2") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 106 105 107 15 sequenceUtil.cpp:354 (set (reg:SF 98)
        (div:SF (reg:SF 99)
            (reg/v:SF 62 [ prob.728 ]))) -1 (nil))

(insn 107 106 108 15 sequenceUtil.cpp:354 (set (reg:DF 100)
        (float_extend:DF (reg:SF 98))) -1 (nil))

(insn 108 107 109 15 sequenceUtil.cpp:354 (set (reg:DF 21 xmm0)
        (reg:DF 100)) -1 (nil))

(call_insn 109 108 110 15 sequenceUtil.cpp:354 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("log") [flags 0x41]  <function_decl 0x7f9013e35c00 log>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
        (nil)))

(insn 110 109 111 15 sequenceUtil.cpp:354 (set (reg:DF 61 [ temp.733 ])
        (reg:DF 21 xmm0)) -1 (nil))

(debug_insn 111 110 112 15 sequenceUtil.cpp:354 (var_location:SF val (float_truncate:SF (div:DF (mult:DF (float_extend:DF (reg/v:SF 62 [ prob.728 ]))
                (reg:DF 61 [ temp.733 ]))
            (const_double:DF 6.9314718055994528622676398299518041312694549560546875e-1 [0x0.b17217f7d1cf78p+0])))) -1 (nil))

(insn 112 111 113 15 sequenceUtil.cpp:355 (set (reg:DF 101)
        (float_extend:DF (reg/v:SF 62 [ prob.728 ]))) -1 (nil))

(insn 113 112 114 15 sequenceUtil.cpp:355 (set (reg:DF 102)
        (mult:DF (reg:DF 101)
            (reg:DF 61 [ temp.733 ]))) -1 (nil))

(insn 114 113 115 15 sequenceUtil.cpp:355 (set (reg:DF 104)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC3") [flags 0x2]) [96 S8 A64])) -1 (expr_list:REG_EQUAL (const_double:DF 6.9314718055994528622676398299518041312694549560546875e-1 [0x0.b17217f7d1cf78p+0])
        (nil)))

(insn 115 114 116 15 sequenceUtil.cpp:355 (set (reg:DF 103)
        (div:DF (reg:DF 102)
            (reg:DF 104))) -1 (nil))

(insn 116 115 117 15 sequenceUtil.cpp:355 (set (reg:SF 105)
        (float_truncate:SF (reg:DF 103))) -1 (nil))

(insn 117 116 118 15 sequenceUtil.cpp:355 (set (reg/v:SF 69 [ entropy ])
        (plus:SF (reg/v:SF 69 [ entropy ])
            (reg:SF 105))) -1 (nil))

(debug_insn 118 117 119 15 sequenceUtil.cpp:355 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))
;; End of basic block 15 -> ( 16)

;; Succ edge  16 [100.0%]  (fallthru)

;; Start of basic block ( 15 14) -> 16
;; Pred edge  15 [100.0%]  (fallthru)
;; Pred edge  14 [100.0%] 
(code_label 119 118 120 16 50 "" [1 uses])

(note 120 119 121 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(debug_insn 121 120 122 16 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))

(debug_insn 122 121 123 16 sequenceUtil.cpp:349 (var_location:SI i (const_int 2 [0x2])) -1 (nil))

(debug_insn 123 122 124 16 (var_location:SI i (const_int 2 [0x2])) -1 (nil))

(debug_insn 124 123 125 16 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))

(insn 125 124 126 16 sequenceUtil.cpp:351 (set (reg:SI 107)
        (sign_extend:SI (mem/s/j:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -14 [0xfffffffffffffff2])) [0 counts+2 S1 A16]))) -1 (nil))

(insn 126 125 127 16 sequenceUtil.cpp:351 (set (reg:SF 106)
        (float:SF (reg:SI 107))) -1 (nil))

(insn 127 126 128 16 sequenceUtil.cpp:351 (set (reg/v:SF 60 [ prob.739 ])
        (div:SF (reg:SF 106)
            (reg:SF 71 [ D.32917 ]))) -1 (nil))

(debug_insn 128 127 129 16 sequenceUtil.cpp:351 (var_location:SF prob (reg/v:SF 60 [ prob.739 ])) -1 (nil))

(insn 129 128 130 16 sequenceUtil.cpp:353 (set (reg:SF 108)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 0.0 [0x0.0p+0])
        (nil)))

(insn 130 129 131 16 sequenceUtil.cpp:353 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:SF 60 [ prob.739 ])
            (reg:SF 108))) -1 (nil))

(jump_insn 131 130 204 16 sequenceUtil.cpp:353 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 211)
            (pc))) 617 {*jcc_1} (expr_list:REG_BR_PROB (const_int 2900 [0xb54])
        (nil)))
;; End of basic block 16 -> ( 18 17)

;; Succ edge  18 [29.0%] 
;; Succ edge  17 [71.0%]  (fallthru)

;; Start of basic block ( 16) -> 17
;; Pred edge  16 [71.0%]  (fallthru)
(note 204 131 132 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(jump_insn 132 204 133 17 sequenceUtil.cpp:353 (set (pc)
        (label_ref 150)) -1 (nil))
;; End of basic block 17 -> ( 19)

;; Succ edge  19 [100.0%] 

(barrier 133 132 211)

;; Start of basic block ( 16) -> 18
;; Pred edge  16 [29.0%] 
(code_label 211 133 135 18 58 "" [1 uses])

(note 135 211 136 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(insn 136 135 137 18 sequenceUtil.cpp:354 (set (reg:SF 110)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC2") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 137 136 138 18 sequenceUtil.cpp:354 (set (reg:SF 109)
        (div:SF (reg:SF 110)
            (reg/v:SF 60 [ prob.739 ]))) -1 (nil))

(insn 138 137 139 18 sequenceUtil.cpp:354 (set (reg:DF 111)
        (float_extend:DF (reg:SF 109))) -1 (nil))

(insn 139 138 140 18 sequenceUtil.cpp:354 (set (reg:DF 21 xmm0)
        (reg:DF 111)) -1 (nil))

(call_insn 140 139 141 18 sequenceUtil.cpp:354 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("log") [flags 0x41]  <function_decl 0x7f9013e35c00 log>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
        (nil)))

(insn 141 140 142 18 sequenceUtil.cpp:354 (set (reg:DF 59 [ temp.744 ])
        (reg:DF 21 xmm0)) -1 (nil))

(debug_insn 142 141 143 18 sequenceUtil.cpp:354 (var_location:SF val (float_truncate:SF (div:DF (mult:DF (float_extend:DF (reg/v:SF 60 [ prob.739 ]))
                (reg:DF 59 [ temp.744 ]))
            (const_double:DF 6.9314718055994528622676398299518041312694549560546875e-1 [0x0.b17217f7d1cf78p+0])))) -1 (nil))

(insn 143 142 144 18 sequenceUtil.cpp:355 (set (reg:DF 112)
        (float_extend:DF (reg/v:SF 60 [ prob.739 ]))) -1 (nil))

(insn 144 143 145 18 sequenceUtil.cpp:355 (set (reg:DF 113)
        (mult:DF (reg:DF 112)
            (reg:DF 59 [ temp.744 ]))) -1 (nil))

(insn 145 144 146 18 sequenceUtil.cpp:355 (set (reg:DF 115)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC3") [flags 0x2]) [96 S8 A64])) -1 (expr_list:REG_EQUAL (const_double:DF 6.9314718055994528622676398299518041312694549560546875e-1 [0x0.b17217f7d1cf78p+0])
        (nil)))

(insn 146 145 147 18 sequenceUtil.cpp:355 (set (reg:DF 114)
        (div:DF (reg:DF 113)
            (reg:DF 115))) -1 (nil))

(insn 147 146 148 18 sequenceUtil.cpp:355 (set (reg:SF 116)
        (float_truncate:SF (reg:DF 114))) -1 (nil))

(insn 148 147 149 18 sequenceUtil.cpp:355 (set (reg/v:SF 69 [ entropy ])
        (plus:SF (reg/v:SF 69 [ entropy ])
            (reg:SF 116))) -1 (nil))

(debug_insn 149 148 150 18 sequenceUtil.cpp:355 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))
;; End of basic block 18 -> ( 19)

;; Succ edge  19 [100.0%]  (fallthru)

;; Start of basic block ( 18 17) -> 19
;; Pred edge  18 [100.0%]  (fallthru)
;; Pred edge  17 [100.0%] 
(code_label 150 149 151 19 52 "" [1 uses])

(note 151 150 152 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(debug_insn 152 151 153 19 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))

(debug_insn 153 152 154 19 sequenceUtil.cpp:349 (var_location:SI i (const_int 3 [0x3])) -1 (nil))

(debug_insn 154 153 155 19 (var_location:SI i (const_int 3 [0x3])) -1 (nil))

(debug_insn 155 154 156 19 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))

(insn 156 155 157 19 sequenceUtil.cpp:351 (set (reg:SI 118)
        (sign_extend:SI (mem/s/j:QI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -13 [0xfffffffffffffff3])) [0 counts+3 S1 A8]))) -1 (nil))

(insn 157 156 158 19 sequenceUtil.cpp:351 (set (reg:SF 117)
        (float:SF (reg:SI 118))) -1 (nil))

(insn 158 157 159 19 sequenceUtil.cpp:351 (set (reg/v:SF 66 [ prob ])
        (div:SF (reg:SF 117)
            (reg:SF 71 [ D.32917 ]))) -1 (nil))

(debug_insn 159 158 160 19 sequenceUtil.cpp:351 (var_location:SF prob (reg/v:SF 66 [ prob ])) -1 (nil))

(insn 160 159 161 19 sequenceUtil.cpp:353 (set (reg:SF 119)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 0.0 [0x0.0p+0])
        (nil)))

(insn 161 160 162 19 sequenceUtil.cpp:353 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:SF 66 [ prob ])
            (reg:SF 119))) -1 (nil))

(jump_insn 162 161 206 19 sequenceUtil.cpp:353 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 212)
            (pc))) 617 {*jcc_1} (expr_list:REG_BR_PROB (const_int 2900 [0xb54])
        (nil)))
;; End of basic block 19 -> ( 21 20)

;; Succ edge  21 [29.0%] 
;; Succ edge  20 [71.0%]  (fallthru)

;; Start of basic block ( 19) -> 20
;; Pred edge  19 [71.0%]  (fallthru)
(note 206 162 163 20 [bb 20] NOTE_INSN_BASIC_BLOCK)

(jump_insn 163 206 164 20 sequenceUtil.cpp:353 (set (pc)
        (label_ref 181)) -1 (nil))
;; End of basic block 20 -> ( 22)

;; Succ edge  22 [100.0%] 

(barrier 164 163 212)

;; Start of basic block ( 19) -> 21
;; Pred edge  19 [29.0%] 
(code_label 212 164 166 21 59 "" [1 uses])

(note 166 212 167 21 [bb 21] NOTE_INSN_BASIC_BLOCK)

(insn 167 166 168 21 sequenceUtil.cpp:354 (set (reg:SF 121)
        (mem/u/c/i:SF (symbol_ref/u:DI ("*.LC2") [flags 0x2]) [95 S4 A32])) -1 (expr_list:REG_EQUAL (const_double:SF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 168 167 169 21 sequenceUtil.cpp:354 (set (reg:SF 120)
        (div:SF (reg:SF 121)
            (reg/v:SF 66 [ prob ]))) -1 (nil))

(insn 169 168 170 21 sequenceUtil.cpp:354 (set (reg:DF 122)
        (float_extend:DF (reg:SF 120))) -1 (nil))

(insn 170 169 171 21 sequenceUtil.cpp:354 (set (reg:DF 21 xmm0)
        (reg:DF 122)) -1 (nil))

(call_insn 171 170 172 21 sequenceUtil.cpp:354 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("log") [flags 0x41]  <function_decl 0x7f9013e35c00 log>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
        (nil)))

(insn 172 171 173 21 sequenceUtil.cpp:354 (set (reg:DF 70 [ D.32924 ])
        (reg:DF 21 xmm0)) -1 (nil))

(debug_insn 173 172 174 21 sequenceUtil.cpp:354 (var_location:SF val (float_truncate:SF (div:DF (mult:DF (float_extend:DF (reg/v:SF 66 [ prob ]))
                (reg:DF 70 [ D.32924 ]))
            (const_double:DF 6.9314718055994528622676398299518041312694549560546875e-1 [0x0.b17217f7d1cf78p+0])))) -1 (nil))

(insn 174 173 175 21 sequenceUtil.cpp:355 (set (reg:DF 123)
        (float_extend:DF (reg/v:SF 66 [ prob ]))) -1 (nil))

(insn 175 174 176 21 sequenceUtil.cpp:355 (set (reg:DF 124)
        (mult:DF (reg:DF 123)
            (reg:DF 70 [ D.32924 ]))) -1 (nil))

(insn 176 175 177 21 sequenceUtil.cpp:355 (set (reg:DF 126)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC3") [flags 0x2]) [96 S8 A64])) -1 (expr_list:REG_EQUAL (const_double:DF 6.9314718055994528622676398299518041312694549560546875e-1 [0x0.b17217f7d1cf78p+0])
        (nil)))

(insn 177 176 178 21 sequenceUtil.cpp:355 (set (reg:DF 125)
        (div:DF (reg:DF 124)
            (reg:DF 126))) -1 (nil))

(insn 178 177 179 21 sequenceUtil.cpp:355 (set (reg:SF 127)
        (float_truncate:SF (reg:DF 125))) -1 (nil))

(insn 179 178 180 21 sequenceUtil.cpp:355 (set (reg/v:SF 69 [ entropy ])
        (plus:SF (reg/v:SF 69 [ entropy ])
            (reg:SF 127))) -1 (nil))

(debug_insn 180 179 181 21 sequenceUtil.cpp:355 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))
;; End of basic block 21 -> ( 22)

;; Succ edge  22 [100.0%]  (fallthru)

;; Start of basic block ( 21 20) -> 22
;; Pred edge  21 [100.0%]  (fallthru)
;; Pred edge  20 [100.0%] 
(code_label 181 180 182 22 54 "" [1 uses])

(note 182 181 183 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(debug_insn 183 182 184 22 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))

(debug_insn 184 183 185 22 sequenceUtil.cpp:349 (var_location:SI i (const_int 4 [0x4])) -1 (nil))

(debug_insn 185 184 186 22 (var_location:SI i (const_int 4 [0x4])) -1 (nil))

(debug_insn 186 185 187 22 (var_location:SF entropy (reg/v:SF 69 [ entropy ])) -1 (nil))

(insn 187 186 191 22 sequenceUtil.cpp:355 (set (reg:SF 72 [ <result> ])
        (reg/v:SF 69 [ entropy ])) -1 (nil))

(insn 191 187 197 22 sequenceUtil.cpp:360 (set (reg/i:SF 21 xmm0)
        (reg:SF 72 [ <result> ])) -1 (nil))

(insn 197 191 0 22 sequenceUtil.cpp:360 (use (reg/i:SF 21 xmm0)) -1 (nil))
;; End of basic block 22 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function std::_Rb_tree_iterator<_Val> std::_Rb_tree<_Key, _Val, _KeyOfValue, _Compare, _Alloc>::_M_insert_(const std::_Rb_tree_node_base*, const std::_Rb_tree_node_base*, const _Val&) [with _Key = char, _Val = std::pair<const char, int>, _KeyOfValue = std::_Select1st<std::pair<const char, int> >, _Compare = std::less<char>, _Alloc = std::allocator<std::pair<const char, int> >] (_ZNSt8_Rb_treeIcSt4pairIKciESt10_Select1stIS2_ESt4lessIcESaIS2_EE10_M_insert_EPKSt18_Rb_tree_node_baseSB_RKS2_)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
Edge 7->9 redirected to 10
merging block 9 into block 8
Merged 8 and 9 without moving.
merging block 11 into block 10
Merged 10 and 11 without moving.
merging block 12 into block 10
Merged 10 and 12 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 7 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 7 1 2 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 2 7 3 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:874 (set (reg/f:DI 64 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil))

(insn 3 2 4 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:874 (set (reg/v/f:DI 65 [ __x ])
        (reg:DI 4 si [ __x ])) -1 (nil))

(insn 4 3 5 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:874 (set (reg/v/f:DI 66 [ __p ])
        (reg:DI 1 dx [ __p ])) -1 (nil))

(insn 5 4 6 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:874 (set (reg/v/f:DI 67 [ __v ])
        (reg:DI 2 cx [ __v ])) -1 (nil))

(note 6 5 9 2 NOTE_INSN_FUNCTION_BEG)

(insn 9 6 10 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 65 [ __x ])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 10 9 11 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 37)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 8500 [0x2134])
        (nil)))
;; End of basic block 2 -> ( 5 3)

;; Succ edge  5 [85.0%] 
;; Succ edge  3 [15.0%]  (fallthru)

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [15.0%]  (fallthru)
(note 11 10 12 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(debug_insn 12 11 13 3 (var_location:DI this (reg/f:DI 64 [ this ])) -1 (nil))

(insn 13 12 14 3 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (parallel [
            (set (reg:DI 68)
                (plus:DI (reg/f:DI 64 [ this ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 14 13 15 3 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 66 [ __p ])
            (reg:DI 68))) -1 (nil))

(jump_insn 15 14 16 3 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 37)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 1500 [0x5dc])
        (nil)))
;; End of basic block 3 -> ( 5 4)

;; Succ edge  5 [15.0%] 
;; Succ edge  4 [85.0%]  (fallthru)

;; Start of basic block ( 3) -> 4
;; Pred edge  3 [85.0%]  (fallthru)
(note 16 15 17 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(debug_insn 17 16 18 4 (var_location:DI __x (reg/v/f:DI 66 [ __p ])) -1 (nil))

(debug_insn 18 17 19 4 (var_location:DI __x (reg/v/f:DI 66 [ __p ])) -1 (nil))

(debug_insn 19 18 20 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:525 (var_location:DI D.4294967178 (reg/v/f:DI 66 [ __p ])) -1 (nil))

(debug_insn 20 19 21 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:529 (var_location:DI D.4294967287 (plus:DI (debug_expr:DI D#118)
        (const_int 32 [0x20]))) -1 (nil))

(debug_insn 21 20 22 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:529 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 22 21 23 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:529 (var_location:DI __x (debug_expr:DI D#9)) -1 (nil))

(debug_insn 23 22 24 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI D.4294967206 (plus:DI (debug_expr:DI D#118)
        (const_int 32 [0x20]))) -1 (nil))

(debug_insn 24 23 25 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 25 24 26 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI __x (reg/v/f:DI 67 [ __v ])) -1 (nil))

(debug_insn 26 25 27 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI D.4294967205 (reg/v/f:DI 67 [ __v ])) -1 (nil))

(debug_insn 27 26 28 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI D.4294967207 (reg/f:DI 64 [ this ])) -1 (nil))

(debug_insn 28 27 29 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI this (debug_expr:DI D#89)) -1 (nil))

(debug_insn 29 28 30 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI __x (debug_expr:DI D#91)) -1 (nil))

(debug_insn 30 29 31 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI __y (debug_expr:DI D#90)) -1 (nil))

(insn 31 30 32 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (reg:QI 69)
        (mem/s:QI (reg/v/f:DI 67 [ __v ]) [0 <variable>.first+0 S1 A32])) -1 (nil))

(insn 32 31 33 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg:QI 69)
            (mem/s:QI (plus:DI (reg/v/f:DI 66 [ __p ])
                    (const_int 32 [0x20])) [0 <variable>._M_value_field.first+0 S1 A64]))) -1 (nil))

(insn 33 32 34 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (reg:QI 70)
        (lt:QI (reg:CCGC 17 flags)
            (const_int 0 [0x0]))) -1 (expr_list:REG_EQUAL (lt:QI (reg:QI 69)
            (mem/s:QI (plus:DI (reg/v/f:DI 66 [ __p ])
                    (const_int 32 [0x20])) [0 <variable>._M_value_field.first+0 S1 A64]))
        (nil)))

(insn 34 33 35 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (parallel [
            (set (reg:SI 58 [ prephitmp.779 ])
                (zero_extend:SI (reg:QI 70)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(jump_insn 35 34 36 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (pc)
        (label_ref 40)) -1 (nil))
;; End of basic block 4 -> ( 6)

;; Succ edge  6 [100.0%] 

(barrier 36 35 37)

;; Start of basic block ( 3 2) -> 5
;; Pred edge  3 [15.0%] 
;; Pred edge  2 [85.0%] 
(code_label 37 36 38 5 67 "" [2 uses])

(note 38 37 39 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 39 38 40 5 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (set (reg:SI 58 [ prephitmp.779 ])
        (const_int 1 [0x1])) -1 (nil))
;; End of basic block 5 -> ( 6)

;; Succ edge  6 [100.0%]  (fallthru)

;; Start of basic block ( 5 4) -> 6
;; Pred edge  5 [100.0%]  (fallthru)
;; Pred edge  4 [100.0%] 
(code_label 40 39 41 6 68 "" [1 uses])

(note 41 40 42 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(debug_insn 42 41 43 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:QI __insert_left (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 43 42 44 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI this (reg/f:DI 64 [ this ])) -1 (nil))

(debug_insn 44 43 45 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:879 (var_location:DI __x (reg/v/f:DI 67 [ __v ])) -1 (nil))

(debug_insn 45 44 46 6 (var_location:DI this (reg/f:DI 64 [ this ])) -1 (nil))

(debug_insn 46 45 47 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:359 (var_location:DI D.4294967209 (reg/f:DI 64 [ this ])) -1 (nil))

(debug_insn 47 46 48 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:359 (var_location:DI D.4294967210 (debug_expr:DI D#87)) -1 (nil))

(debug_insn 48 47 49 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:359 (var_location:DI this (debug_expr:DI D#86)) -1 (nil))

(debug_insn 49 48 50 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:359 (var_location:DI __n (const_int 1 [0x1])) -1 (nil))

(debug_insn 50 49 51 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:359 (var_location:DI D.40531 (const_int 0 [0x0])) -1 (nil))

(debug_insn 51 50 52 6 (var_location:DI this (debug_expr:DI D#86)) -1 (nil))

(insn 52 51 53 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:89 (set (reg:DI 5 di)
        (const_int 40 [0x28])) -1 (nil))

(call_insn 53 52 54 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:89 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_Znwm") [flags 0x41]  <function_decl 0x7f9013d1e500 operator new>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 54 53 55 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:89 (set (reg/f:DI 60 [ D.40538 ])
        (reg:DI 0 ax)) -1 (nil))

(insn 55 54 56 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:89 (set (reg/v/f:DI 61 [ __z ])
        (reg/f:DI 60 [ D.40538 ])) -1 (nil))

(debug_insn 56 55 57 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:369 (var_location:DI __tmp (reg/v/f:DI 61 [ __z ])) -1 (nil))

(debug_insn 57 56 58 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:371 (var_location:DI D.4294967208 (plus:DI (reg/v/f:DI 61 [ __z ])
        (const_int 32 [0x20]))) -1 (nil))

(debug_insn 58 57 59 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:371 (var_location:DI this (reg/f:DI 64 [ this ])) -1 (nil))

(debug_insn 59 58 60 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:354 (var_location:DI this (reg/f:DI 64 [ this ])) -1 (nil))

(debug_insn 60 59 61 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:350 (var_location:DI D.4294967271 (reg/f:DI 64 [ this ])) -1 (nil))

(debug_insn 61 60 62 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:354 (var_location:DI this (debug_expr:DI D#26)) -1 (nil))

(debug_insn 62 61 63 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:354 (var_location:DI D.39291 (debug_expr:DI D#25)) -1 (nil))

(debug_insn 63 62 64 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/allocator.h:107 (var_location:DI this (debug_expr:DI D#26)) -1 (nil))

(debug_insn 64 63 65 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:371 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 65 64 66 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:371 (var_location:DI __p (debug_expr:DI D#88)) -1 (nil))

(debug_insn 66 65 67 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:371 (var_location:DI __val (reg/v/f:DI 67 [ __v ])) -1 (nil))

(debug_insn 67 66 68 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:105 (var_location:DI D.4294967177 (plus:DI (reg/v/f:DI 61 [ __z ])
        (const_int 32 [0x20]))) -1 (nil))

(debug_insn 68 67 69 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:105 (var_location:DI D.40478 (const_int 8 [0x8])) -1 (nil))

(debug_insn 69 68 70 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:105 (var_location:DI __p (debug_expr:DI D#119)) -1 (nil))

(insn 70 69 71 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:105 (parallel [
            (set (reg:DI 71)
                (plus:DI (reg/v/f:DI 61 [ __z ])
                    (const_int 32 [0x20])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 71 70 72 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:105 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 71)
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 72 71 73 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:105 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 78)
            (pc))) 617 {*jcc_1} (expr_list:REG_BR_PROB (const_int 1014 [0x3f6])
        (nil)))
;; End of basic block 6 -> ( 7 8)

;; Succ edge  7 [89.9%]  (fallthru)
;; Succ edge  8 [10.1%] 

;; Start of basic block ( 6) -> 7
;; Pred edge  6 [89.9%]  (fallthru)
(note 73 72 74 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 74 73 75 7 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:105 (set (reg:DI 72)
        (mem/s:DI (reg/v/f:DI 67 [ __v ]) [72 S8 A32])) -1 (nil))

(insn 75 74 78 7 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/ext/new_allocator.h:105 (set (mem/s:DI (plus:DI (reg/v/f:DI 61 [ __z ])
                (const_int 32 [0x20])) [72 <variable>._M_value_field+0 S8 A64])
        (reg:DI 72)) -1 (nil))
;; End of basic block 7 -> ( 8)

;; Succ edge  8 [100.0%]  (fallthru)

;; Start of basic block ( 7 6) -> 8
;; Pred edge  7 [100.0%]  (fallthru)
;; Pred edge  6 [10.1%] 
(code_label 78 75 105 8 70 "" [1 uses])

(note 105 78 79 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(debug_insn 79 105 80 8 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 80 79 81 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/allocator.h:109 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 81 80 82 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:881 (var_location:DI __z (reg/v/f:DI 61 [ __z ])) -1 (nil))

(insn 82 81 83 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:883 (set (reg/f:DI 59 [ SR.765 ])
        (reg/v/f:DI 61 [ __z ])) -1 (nil))

(insn 83 82 84 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:883 (parallel [
            (set (reg:DI 73)
                (plus:DI (reg/f:DI 64 [ this ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 84 83 85 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:883 (set (reg:DI 2 cx)
        (reg:DI 73)) -1 (nil))

(insn 85 84 86 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:883 (set (reg:DI 1 dx)
        (reg/v/f:DI 66 [ __p ])) -1 (nil))

(insn 86 85 87 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:883 (set (reg:DI 4 si)
        (reg/f:DI 59 [ SR.765 ])) -1 (nil))

(insn 87 86 88 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:883 (set (reg:SI 5 di)
        (reg:SI 58 [ prephitmp.779 ])) -1 (nil))

(call_insn 88 87 89 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:883 (call (mem:QI (symbol_ref:DI ("_ZSt29_Rb_tree_insert_and_rebalancebPSt18_Rb_tree_node_baseS0_RS_") [flags 0x41]  <function_decl 0x7f90127a9000 _Rb_tree_insert_and_rebalance>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))

(insn 89 88 90 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:886 (parallel [
            (set (mem/s:DI (plus:DI (reg/f:DI 64 [ this ])
                        (const_int 40 [0x28])) [14 <variable>._M_impl._M_node_count+0 S8 A64])
                (plus:DI (mem/s:DI (plus:DI (reg/f:DI 64 [ this ])
                            (const_int 40 [0x28])) [14 <variable>._M_impl._M_node_count+0 S8 A64])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 90 89 91 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:886 (var_location:DI this (clobber (const_int 0 [0x0]))) -1 (nil))

(debug_insn 91 90 92 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:886 (var_location:DI __x (reg/v/f:DI 61 [ __z ])) -1 (nil))

(insn 92 91 93 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:886 (set (reg:DI 62 [ D.37494 ])
        (reg/f:DI 59 [ SR.765 ])) -1 (nil))

(insn 93 92 97 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:886 (set (reg:DI 63 [ <result> ])
        (reg:DI 62 [ D.37494 ])) -1 (nil))

(insn 97 93 103 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:888 (set (reg/i:DI 0 ax)
        (reg:DI 63 [ <result> ])) -1 (nil))

(insn 103 97 0 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/stl_tree.h:888 (use (reg/i:DI 0 ax)) -1 (nil))
;; End of basic block 8 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function std::basic_string<_CharT, _Traits, _Alloc> std::operator+(const std::basic_string<_CharT, _Traits, _Alloc>&, const _CharT*) [with _CharT = char, _Traits = std::char_traits<char>, _Alloc = std::allocator<char>] (_ZStplIcSt11char_traitsIcESaIcEESbIT_T0_T1_ERKS6_PKS3_)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
Forwarding edge 2->4 to 6 failed.
Forwarding edge 2->4 to 6 failed.
merging block 7 into block 6
Merged 6 and 7 without moving.


try_optimize_cfg iteration 2

Forwarding edge 2->4 to 6 failed.


try_optimize_cfg iteration 1

Forwarding edge 2->3 to 6 failed.
merging block 5 into block 4
Merged 4 and 5 without moving.


try_optimize_cfg iteration 2

Forwarding edge 2->3 to 6 failed.
(note 1 0 6 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 6 1 2 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 2 6 3 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2198 (set (reg/f:DI 60 [ <result> ])
        (reg:DI 5 di [ D.42204 ])) -1 (nil))

(insn 3 2 4 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2198 (set (reg/v/f:DI 61 [ __lhs ])
        (reg:DI 4 si [ __lhs ])) -1 (nil))

(insn 4 3 5 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2198 (set (reg/v/f:DI 62 [ __rhs ])
        (reg:DI 1 dx [ __rhs ])) -1 (nil))

(note 5 4 8 2 NOTE_INSN_FUNCTION_BEG)

(insn 8 5 9 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2201 (set (reg/f:DI 59 [ __str.131 ])
        (reg/f:DI 60 [ <result> ])) -1 (nil))

(insn 9 8 10 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2201 (set (reg:DI 4 si)
        (reg/v/f:DI 61 [ __lhs ])) -1 (nil))

(insn 10 9 11 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2201 (set (reg:DI 5 di)
        (reg/f:DI 59 [ __str.131 ])) -1 (nil))

(call_insn 11 10 12 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2201 (call (mem:QI (symbol_ref:DI ("_ZNSsC1ERKSs") [flags 0x41]  <function_decl 0x7f901329b400 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(debug_insn 12 11 13 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2202 (var_location:DI this (reg/f:DI 59 [ __str.131 ])) -1 (nil))

(debug_insn 13 12 14 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2202 (var_location:DI __s (reg/v/f:DI 62 [ __rhs ])) -1 (nil))

(debug_insn 14 13 15 2 (var_location:DI __s (reg/v/f:DI 62 [ __rhs ])) -1 (nil))

(insn 15 14 16 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:263 (set (reg:DI 5 di)
        (reg/v/f:DI 62 [ __rhs ])) -1 (nil))

(call_insn/i 16 15 17 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:263 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strlen") [flags 0x41]  <function_decl 0x7f9013e83000 __builtin_strlen>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 17 16 18 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:263 (set (reg:DI 64)
        (reg:DI 0 ax)) -1 (nil))

(insn 18 17 19 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:263 (set (reg:DI 65)
        (reg:DI 64)) -1 (nil))

(insn 19 18 20 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/char_traits.h:263 (set (reg:DI 58 [ D.42178 ])
        (reg:DI 65)) -1 (nil))

(insn 20 19 21 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:871 (set (reg:DI 1 dx)
        (reg:DI 58 [ D.42178 ])) -1 (nil))

(insn 21 20 22 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:871 (set (reg:DI 4 si)
        (reg/v/f:DI 62 [ __rhs ])) -1 (nil))

(insn 22 21 23 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:871 (set (reg:DI 5 di)
        (reg/f:DI 59 [ __str.131 ])) -1 (nil))

(call_insn 23 22 24 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:871 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6appendEPKcm") [flags 0x41]  <function_decl 0x7f901324a600 append>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 2 -> ( 3 4)

;; Succ edge  3 [100.0%]  (fallthru)
;; Succ edge  4 (ab,abcall,eh)

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [100.0%]  (fallthru)
(note 24 23 25 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(jump_insn 25 24 26 3 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:871 (set (pc)
        (label_ref 33)) -1 (nil))
;; End of basic block 3 -> ( 6)

;; Succ edge  6 [100.0%] 

(barrier 26 25 48)

;; Start of basic block ( 2) -> 4
;; Pred edge  2 (ab,abcall,eh)
(code_label/s 48 26 51 4 77 "" [1 uses])

(note 51 48 49 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 49 51 50 4 (set (reg:DI 66)
        (reg:DI 0 ax)) -1 (nil))

(insn 50 49 27 4 (set (reg:DI 67)
        (reg:DI 1 dx)) -1 (nil))

(note/s 27 50 29 4 "" NOTE_INSN_DELETED_LABEL 75)

(insn 29 27 30 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2203 (set (reg:DI 5 di)
        (reg/f:DI 59 [ __str.131 ])) -1 (nil))

(call_insn 30 29 45 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2203 (call (mem:QI (symbol_ref:DI ("_ZNSsD1Ev") [flags 0x41]  <function_decl 0x7f901329ba00 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (expr_list:REG_EH_REGION (const_int 2 [0x2])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 45 30 46 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2203 (set (reg:DI 5 di)
        (reg:DI 66)) -1 (nil))

(call_insn 46 45 32 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2203 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 4 -> ()


(barrier 32 46 33)

;; Start of basic block ( 3) -> 6
;; Pred edge  3 [100.0%] 
(code_label 33 32 43 6 74 "" [1 uses])

(note 43 33 34 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 34 43 35 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2204 (set (reg/i:DI 0 ax)
        (reg/f:DI 60 [ <result> ])) -1 (nil))

(insn 35 34 41 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2204 (set (reg/i:DI 0 ax)
        (reg/f:DI 60 [ <result> ])) -1 (nil))

(insn 41 35 0 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:2204 (use (reg/i:DI 0 ax)) -1 (nil))
;; End of basic block 6 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function std::string decode_kmer_from_intval(kmer_int_type_t, unsigned int) (_Z23decode_kmer_from_intvalyj)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
Forwarding edge 4->5 to 17 failed.
Redirecting jump 55 from 16 to 17.
merging block 7 into block 6
Merged 6 and 7 without moving.
merging block 8 into block 6
Merged 6 and 8 without moving.
Edge 13->16 redirected to 17
deleting block 16
merging block 18 into block 17
Merged 17 and 18 without moving.


try_optimize_cfg iteration 2

Forwarding edge 4->5 to 17 failed.


try_optimize_cfg iteration 1

Forwarding edge 3->4 to 15 failed.
merging block 6 into block 5
Merged 5 and 6 without moving.
merging block 14 into block 13
Merged 13 and 14 without moving.


try_optimize_cfg iteration 2

Forwarding edge 3->4 to 15 failed.
(note 33 0 38 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 38 33 34 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 34 38 35 2 sequenceUtil.cpp:304 (set (reg/f:DI 95 [ <result> ])
        (reg:DI 5 di [ D.42262 ])) -1 (nil))

(insn 35 34 36 2 sequenceUtil.cpp:304 (set (reg/v:DI 96 [ intval ])
        (reg:DI 4 si [ intval ])) -1 (nil))

(insn 36 35 37 2 sequenceUtil.cpp:304 (set (reg/v:SI 97 [ kmer_length ])
        (reg:SI 1 dx [ kmer_length ])) -1 (nil))

(note 37 36 40 2 NOTE_INSN_FUNCTION_BEG)

(debug_insn 40 37 41 2 (var_location:DI this (plus:DI (reg/f:DI 54 virtual-stack-vars)
        (const_int -1 [0xffffffffffffffff]))) -1 (nil))

(debug_insn 41 40 42 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/allocator.h:101 (var_location:DI this (plus:DI (reg/f:DI 54 virtual-stack-vars)
        (const_int -1 [0xffffffffffffffff]))) -1 (nil))

(insn 42 41 43 2 sequenceUtil.cpp:306 (set (reg/f:DI 94 [ kmer.100 ])
        (reg/f:DI 95 [ <result> ])) -1 (nil))

(insn 43 42 44 2 sequenceUtil.cpp:306 (parallel [
            (set (reg:DI 98)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 44 43 45 2 sequenceUtil.cpp:306 (set (reg:DI 99)
        (zero_extend:DI (reg/v:SI 97 [ kmer_length ]))) -1 (nil))

(insn 45 44 46 2 sequenceUtil.cpp:306 (set (reg:DI 2 cx)
        (reg:DI 98)) -1 (nil))

(insn 46 45 47 2 sequenceUtil.cpp:306 (set (reg:SI 1 dx)
        (const_int 32 [0x20])) -1 (nil))

(insn 47 46 48 2 sequenceUtil.cpp:306 (set (reg:DI 4 si)
        (reg:DI 99)) -1 (nil))

(insn 48 47 49 2 sequenceUtil.cpp:306 (set (reg:DI 5 di)
        (reg/f:DI 94 [ kmer.100 ])) -1 (nil))

(call_insn 49 48 50 2 sequenceUtil.cpp:306 (call (mem:QI (symbol_ref:DI ("_ZNSsC1EmcRKSaIcE") [flags 0x41]  <function_decl 0x7f9013291a00 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:SI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 2 -> ( 3 5)

;; Succ edge  3 [100.0%]  (fallthru)
;; Succ edge  5 (ab,abcall,eh)

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [100.0%]  (fallthru)
(note 50 49 51 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(debug_insn 51 50 52 3 (var_location:SI i (const_int 1 [0x1])) -1 (nil))

(debug_insn 52 51 53 3 (var_location:DI intval (reg/v:DI 96 [ intval ])) -1 (nil))

(insn 53 52 54 3 sequenceUtil.cpp:308 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 97 [ kmer_length ])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 54 53 57 3 sequenceUtil.cpp:308 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9550 [0x254e])
        (nil)))
;; End of basic block 3 -> ( 7 4)

;; Succ edge  7 [95.5%] 
;; Succ edge  4 [4.5%]  (fallthru)

;; Start of basic block ( 3) -> 4
;; Pred edge  3 [4.5%]  (fallthru)
(note 57 54 55 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(jump_insn 55 57 56 4 sequenceUtil.cpp:308 (set (pc)
        (label_ref:DI 125)) 638 {jump} (nil))
;; End of basic block 4 -> ( 15)

;; Succ edge  15 [100.0%] 

(barrier 56 55 148)

;; Start of basic block ( 2) -> 5
;; Pred edge  2 (ab,abcall,eh)
(code_label/s 148 56 151 5 90 "" [1 uses])

(note 151 148 149 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 149 151 150 5 (set (reg:DI 101)
        (reg:DI 0 ax)) -1 (nil))

(insn 150 149 58 5 (set (reg:DI 100)
        (reg:DI 1 dx)) -1 (nil))

(note/s 58 150 60 5 "" NOTE_INSN_DELETED_LABEL 82)

(insn 60 58 61 5 sequenceUtil.cpp:308 (set (reg:SI 91 [ save_filt.371 ])
        (subreg:SI (reg:DI 100) 0)) -1 (nil))

(insn 61 60 62 5 sequenceUtil.cpp:308 (set (reg/f:DI 92 [ save_eptr.370 ])
        (reg:DI 101)) -1 (nil))

(debug_insn 62 61 63 5 (var_location:DI this (plus:DI (reg/f:DI 54 virtual-stack-vars)
        (const_int -1 [0xffffffffffffffff]))) -1 (nil))

(debug_insn 63 62 67 5 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/allocator.h:109 (var_location:DI this (plus:DI (reg/f:DI 54 virtual-stack-vars)
        (const_int -1 [0xffffffffffffffff]))) -1 (nil))

(insn 67 63 68 5 sequenceUtil.cpp:308 (set (reg:DI 101)
        (reg/f:DI 92 [ save_eptr.370 ])) -1 (nil))

(insn 68 67 141 5 sequenceUtil.cpp:308 (set (reg:DI 100)
        (sign_extend:DI (reg:SI 91 [ save_filt.371 ]))) -1 (nil))

(insn 141 68 142 5 sequenceUtil.cpp:308 (set (reg:DI 5 di)
        (reg:DI 101)) -1 (nil))

(call_insn 142 141 70 5 sequenceUtil.cpp:308 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 5 -> ()


(barrier 70 142 71)

;; Start of basic block ( 3) -> 7
;; Pred edge  3 [95.5%] 
(code_label 71 70 72 7 80 "" [1 uses])

(note 72 71 73 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 73 72 112 7 sequenceUtil.cpp:308 (set (reg/v:SI 93 [ i ])
        (const_int 1 [0x1])) -1 (nil))
;; End of basic block 7 -> ( 8)

;; Succ edge  8 [100.0%]  (fallthru)

;; Start of basic block ( 7 12) -> 8
;; Pred edge  7 [100.0%]  (fallthru)
;; Pred edge  12 [100.0%]  (dfs_back)
(code_label 112 73 74 8 86 "" [1 uses])

(note 74 112 75 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(debug_insn 75 74 76 8 sequenceUtil.cpp:310 (var_location:SI D.4294967174 (subreg:SI (reg/v:DI 96 [ intval ]) 0)) -1 (nil))

(debug_insn 76 75 77 8 sequenceUtil.cpp:310 (var_location:SI base_num (and:SI (debug_expr:SI D#122)
        (const_int 3 [0x3]))) -1 (nil))

(debug_insn 77 76 78 8 sequenceUtil.cpp:312 (var_location:SI D.4294967175 (minus:SI (reg/v:SI 97 [ kmer_length ])
        (reg/v:SI 93 [ i ]))) -1 (nil))

(debug_insn 78 77 79 8 sequenceUtil.cpp:312 (var_location:DI this (reg/f:DI 94 [ kmer.100 ])) -1 (nil))

(debug_insn 79 78 80 8 sequenceUtil.cpp:312 (var_location:DI __pos (zero_extend:DI (debug_expr:SI D#121))) -1 (nil))

(debug_insn 80 79 81 8 sequenceUtil.cpp:312 (var_location:DI this (reg/f:DI 94 [ kmer.100 ])) -1 (nil))

(debug_insn 81 80 82 8 (var_location:DI this (reg/f:DI 94 [ kmer.100 ])) -1 (nil))

(debug_insn 82 81 83 8 (var_location:DI this (reg/f:DI 94 [ kmer.100 ])) -1 (nil))

(insn 83 82 84 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:278 (set (reg/f:DI 90 [ prephitmp.842 ])
        (mem/s/f:DI (reg/f:DI 94 [ kmer.100 ]) [32 <variable>._M_dataplus._M_p+0 S8 A64])) -1 (nil))

(debug_insn 84 83 85 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:301 (var_location:DI this (plus:DI (reg/f:DI 90 [ prephitmp.842 ])
        (const_int -24 [0xffffffffffffffe8]))) -1 (nil))

(insn 85 84 86 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:301 (set (reg:CCGOC 17 flags)
        (compare:CCGOC (mem/s:SI (plus:DI (reg/f:DI 90 [ prephitmp.842 ])
                    (const_int -8 [0xfffffffffffffff8])) [5 <variable>.D.11486._M_refcount+0 S4 A64])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 86 85 87 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:301 (set (pc)
        (if_then_else (lt (reg:CCGOC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 92)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 3666 [0xe52])
        (nil)))
;; End of basic block 8 -> ( 9 11)

;; Succ edge  9 [63.3%]  (fallthru)
;; Succ edge  11 [36.7%] 

;; Start of basic block ( 8) -> 9
;; Pred edge  8 [63.3%]  (fallthru)
(note 87 86 88 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 88 87 89 9 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:302 (set (reg:DI 5 di)
        (reg/f:DI 94 [ kmer.100 ])) -1 (nil))

(call_insn 89 88 90 9 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:302 (call (mem:QI (symbol_ref:DI ("_ZNSs12_M_leak_hardEv") [flags 0x41]  <function_decl 0x7f901322ec00 _M_leak_hard>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (expr_list:REG_EH_REGION (const_int 3 [0x3])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 9 -> ( 10 13)

;; Succ edge  10 [100.0%]  (fallthru)
;; Succ edge  13 (ab,abcall,eh)

;; Start of basic block ( 9) -> 10
;; Pred edge  9 [100.0%]  (fallthru)
(note 90 89 91 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(insn 91 90 92 10 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:302 (set (reg/f:DI 90 [ prephitmp.842 ])
        (mem/s/f:DI (reg/f:DI 94 [ kmer.100 ]) [32 <variable>._M_dataplus._M_p+0 S8 A64])) -1 (nil))
;; End of basic block 10 -> ( 11)

;; Succ edge  11 [100.0%]  (fallthru)

;; Start of basic block ( 8 10) -> 11
;; Pred edge  8 [36.7%] 
;; Pred edge  10 [100.0%]  (fallthru)
(code_label 92 91 93 11 85 "" [1 uses])

(note 93 92 94 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(debug_insn 94 93 95 11 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:746 (var_location:DI this (reg/f:DI 94 [ kmer.100 ])) -1 (nil))

(insn 95 94 96 11 sequenceUtil.cpp:312 (parallel [
            (set (reg:SI 102)
                (minus:SI (reg/v:SI 97 [ kmer_length ])
                    (reg/v:SI 93 [ i ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 96 95 97 11 sequenceUtil.cpp:312 (set (reg:DI 103)
        (zero_extend:DI (reg:SI 102))) -1 (nil))

(insn 97 96 98 11 sequenceUtil.cpp:312 (set (reg:DI 104)
        (sign_extend:DI (subreg:SI (reg/v:DI 96 [ intval ]) 0))) -1 (nil))

(insn 98 97 99 11 sequenceUtil.cpp:312 (parallel [
            (set (reg:DI 105)
                (and:DI (reg:DI 104)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 99 98 100 11 sequenceUtil.cpp:312 (set (reg/f:DI 106)
        (symbol_ref:DI ("_int_to_base") [flags 0x2]  <var_decl 0x7f9012470b40 _int_to_base>)) -1 (nil))

(insn 100 99 101 11 sequenceUtil.cpp:312 (set (reg:QI 107)
        (mem/s/j:QI (plus:DI (reg/f:DI 106)
                (reg:DI 105)) [0 _int_to_base S1 A8])) -1 (nil))

(insn 101 100 102 11 sequenceUtil.cpp:312 (set (mem:QI (plus:DI (reg/f:DI 90 [ prephitmp.842 ])
                (reg:DI 103)) [0 S1 A8])
        (reg:QI 107)) -1 (nil))

(debug_insn 102 101 103 11 sequenceUtil.cpp:316 (var_location:DI D.4294967176 (lshiftrt:DI (reg/v:DI 96 [ intval ])
        (const_int 2 [0x2]))) -1 (nil))

(debug_insn 103 102 104 11 sequenceUtil.cpp:316 (var_location:DI intval (debug_expr:DI D#120)) -1 (nil))

(insn 104 103 105 11 sequenceUtil.cpp:308 (parallel [
            (set (reg/v:SI 93 [ i ])
                (plus:SI (reg/v:SI 93 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 105 104 106 11 sequenceUtil.cpp:308 (var_location:SI i (reg/v:SI 93 [ i ])) -1 (nil))

(debug_insn 106 105 107 11 (var_location:SI i (reg/v:SI 93 [ i ])) -1 (nil))

(debug_insn 107 106 108 11 (var_location:DI intval (debug_expr:DI D#120)) -1 (nil))

(insn 108 107 109 11 sequenceUtil.cpp:308 (set (reg:CC 17 flags)
        (compare:CC (reg/v:SI 97 [ kmer_length ])
            (reg/v:SI 93 [ i ]))) -1 (nil))

(jump_insn 109 108 110 11 sequenceUtil.cpp:308 (set (pc)
        (if_then_else (ltu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 125)
            (pc))) 617 {*jcc_1} (expr_list:REG_BR_PROB (const_int 450 [0x1c2])
        (nil)))
;; End of basic block 11 -> ( 12 15)

;; Succ edge  12 [95.5%]  (fallthru)
;; Succ edge  15 [4.5%] 

;; Start of basic block ( 11) -> 12
;; Pred edge  11 [95.5%]  (fallthru)
(note 110 109 111 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(insn 111 110 113 12 sequenceUtil.cpp:316 (parallel [
            (set (reg/v:DI 96 [ intval ])
                (lshiftrt:DI (reg/v:DI 96 [ intval ])
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(jump_insn 113 111 114 12 sequenceUtil.cpp:316 (set (pc)
        (label_ref 112)) -1 (nil))
;; End of basic block 12 -> ( 8)

;; Succ edge  8 [100.0%]  (dfs_back)

(barrier 114 113 144)

;; Start of basic block ( 9) -> 13
;; Pred edge  9 (ab,abcall,eh)
(code_label/s 144 114 147 13 89 "" [1 uses])

(note 147 144 145 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(insn 145 147 146 13 (set (reg:DI 101)
        (reg:DI 0 ax)) -1 (nil))

(insn 146 145 115 13 (set (reg:DI 100)
        (reg:DI 1 dx)) -1 (nil))

(note/s 115 146 117 13 "" NOTE_INSN_DELETED_LABEL 87)

(insn 117 115 118 13 sequenceUtil.cpp:319 (set (reg:DI 5 di)
        (reg/f:DI 94 [ kmer.100 ])) -1 (nil))

(call_insn 118 117 138 13 sequenceUtil.cpp:319 (call (mem:QI (symbol_ref:DI ("_ZNSsD1Ev") [flags 0x41]  <function_decl 0x7f901329ba00 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (expr_list:REG_EH_REGION (const_int 4 [0x4])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 138 118 139 13 sequenceUtil.cpp:319 (set (reg:DI 5 di)
        (reg:DI 101)) -1 (nil))

(call_insn 139 138 120 13 sequenceUtil.cpp:319 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 13 -> ()


(barrier 120 139 125)

;; Start of basic block ( 11 4) -> 15
;; Pred edge  11 [4.5%] 
;; Pred edge  4 [100.0%] 
(code_label 125 120 136 15 79 "" [2 uses])

(note 136 125 126 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(insn 126 136 127 15 sequenceUtil.cpp:320 (set (reg/i:DI 0 ax)
        (reg/f:DI 95 [ <result> ])) -1 (nil))

(insn 127 126 133 15 sequenceUtil.cpp:320 (set (reg/i:DI 0 ax)
        (reg/f:DI 95 [ <result> ])) -1 (nil))

(insn 133 127 0 15 sequenceUtil.cpp:320 (use (reg/i:DI 0 ax)) -1 (nil))
;; End of basic block 15 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)


;; Function bool contains_non_gatc(std::string) (_Z17contains_non_gatcSs)



try_optimize_cfg iteration 1

merging block 3 into block 2
Merged 2 and 3 without moving.
merging block 12 into block 11
Merged 11 and 12 without moving.
merging block 13 into block 11
Merged 11 and 13 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 2 0 5 NOTE_INSN_DELETED)

;; Start of basic block ( 0) -> 2
;; Pred edge  ENTRY [100.0%]  (fallthru)
(note 5 2 3 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 3 5 4 2 sequenceUtil.cpp:30 (set (reg/v/f:DI 65 [ kmer ])
        (reg:DI 5 di [ kmer ])) -1 (nil))

(note 4 3 7 2 NOTE_INSN_FUNCTION_BEG)

(debug_insn 7 4 8 2 sequenceUtil.cpp:32 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 8 7 9 2 (var_location:SI i (const_int 0 [0x0])) -1 (nil))

(debug_insn 9 8 10 2 sequenceUtil.cpp:32 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 10 9 11 2 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 11 10 12 2 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(insn 12 11 13 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:278 (set (reg/f:DI 59 [ prephitmp.887 ])
        (mem/s/f:DI (reg/v/f:DI 65 [ kmer ]) [32 <variable>._M_dataplus._M_p+0 S8 A64])) -1 (nil))

(insn 13 12 14 2 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:286 (parallel [
            (set (reg/f:DI 60 [ D.39717 ])
                (plus:DI (reg/f:DI 59 [ prephitmp.887 ])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 14 13 15 2 sequenceUtil.cpp:32 (set (reg:CCZ 17 flags)
        (compare:CCZ (mem/s:DI (reg/f:DI 60 [ D.39717 ]) [14 <variable>.D.11486._M_length+0 S8 A64])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 15 14 16 2 sequenceUtil.cpp:32 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 58)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 450 [0x1c2])
        (nil)))
;; End of basic block 2 -> ( 3 9)

;; Succ edge  3 [95.5%]  (fallthru)
;; Succ edge  9 [4.5%] 

;; Start of basic block ( 2) -> 3
;; Pred edge  2 [95.5%]  (fallthru)
(note 16 15 17 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 17 16 18 3 sequenceUtil.cpp:32 (set (reg:DI 63 [ D.31789 ])
        (const_int 0 [0x0])) -1 (nil))

(insn 18 17 55 3 sequenceUtil.cpp:32 (set (reg/v:SI 61 [ i ])
        (const_int 0 [0x0])) -1 (nil))
;; End of basic block 3 -> ( 4)

;; Succ edge  4 [100.0%]  (fallthru)

;; Start of basic block ( 3 8) -> 4
;; Pred edge  3 [100.0%]  (fallthru)
;; Pred edge  8 [95.5%]  (dfs_back)
(code_label 55 18 19 4 98 "" [1 uses])

(note 19 55 20 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(debug_insn 20 19 21 4 sequenceUtil.cpp:33 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 21 20 22 4 sequenceUtil.cpp:33 (var_location:DI __pos (reg:DI 63 [ D.31789 ])) -1 (nil))

(debug_insn 22 21 23 4 sequenceUtil.cpp:33 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 23 22 24 4 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 24 23 25 4 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 25 24 26 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:301 (var_location:DI this (reg/f:DI 60 [ D.39717 ])) -1 (nil))

(insn 26 25 27 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:301 (set (reg:CCGOC 17 flags)
        (compare:CCGOC (mem/s:SI (plus:DI (reg/f:DI 60 [ D.39717 ])
                    (const_int 16 [0x10])) [5 <variable>.D.11486._M_refcount+0 S4 A64])
            (const_int 0 [0x0]))) -1 (nil))

(jump_insn 27 26 28 4 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:301 (set (pc)
        (if_then_else (lt (reg:CCGOC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 32)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 3666 [0xe52])
        (nil)))
;; End of basic block 4 -> ( 5 6)

;; Succ edge  5 [63.3%]  (fallthru)
;; Succ edge  6 [36.7%] 

;; Start of basic block ( 4) -> 5
;; Pred edge  4 [63.3%]  (fallthru)
(note 28 27 29 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 29 28 30 5 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:302 (set (reg:DI 5 di)
        (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(call_insn 30 29 31 5 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:302 (call (mem:QI (symbol_ref:DI ("_ZNSs12_M_leak_hardEv") [flags 0x41]  <function_decl 0x7f901322ec00 _M_leak_hard>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 31 30 32 5 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:302 (set (reg/f:DI 59 [ prephitmp.887 ])
        (mem/s/f:DI (reg/v/f:DI 65 [ kmer ]) [32 <variable>._M_dataplus._M_p+0 S8 A64])) -1 (nil))
;; End of basic block 5 -> ( 6)

;; Succ edge  6 [100.0%]  (fallthru)

;; Start of basic block ( 5 4) -> 6
;; Pred edge  5 [100.0%]  (fallthru)
;; Pred edge  4 [36.7%] 
(code_label 32 31 33 6 95 "" [1 uses])

(note 33 32 34 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(debug_insn 34 33 35 6 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:746 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 35 34 36 6 sequenceUtil.cpp:33 (var_location:QI D.4294967173 (mem:QI (plus:DI (reg/f:DI 59 [ prephitmp.887 ])
            (reg:DI 63 [ D.31789 ])) [0 S1 A8])) -1 (nil))

(debug_insn 36 35 37 6 sequenceUtil.cpp:33 (var_location:QI c (debug_expr:QI D#123)) -1 (nil))

(insn 37 36 38 6 sequenceUtil.cpp:35 (set (reg:DI 66)
        (sign_extend:DI (mem:QI (plus:DI (reg/f:DI 59 [ prephitmp.887 ])
                    (reg:DI 63 [ D.31789 ])) [0 S1 A8]))) -1 (nil))

(insn 38 37 39 6 sequenceUtil.cpp:35 (set (reg/f:DI 67)
        (symbol_ref:DI ("_base_to_int") [flags 0x2]  <var_decl 0x7f9012470be0 _base_to_int>)) -1 (nil))

(insn 39 38 40 6 sequenceUtil.cpp:35 (set (reg:CC 17 flags)
        (compare:CC (mem/s/j:QI (plus:DI (reg/f:DI 67)
                    (reg:DI 66)) [0 _base_to_int S1 A8])
            (const_int 3 [0x3]))) -1 (nil))

(jump_insn 40 39 41 6 sequenceUtil.cpp:35 (set (pc)
        (if_then_else (leu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 45)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9550 [0x254e])
        (nil)))
;; End of basic block 6 -> ( 7 8)

;; Succ edge  7 [4.5%]  (fallthru)
;; Succ edge  8 [95.5%] 

;; Start of basic block ( 6) -> 7
;; Pred edge  6 [4.5%]  (fallthru)
(note 41 40 42 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 42 41 43 7 sequenceUtil.cpp:35 (set (reg:QI 62 [ D.31798 ])
        (const_int 1 [0x1])) -1 (nil))

(jump_insn 43 42 44 7 sequenceUtil.cpp:35 (set (pc)
        (label_ref 61)) -1 (nil))
;; End of basic block 7 -> ( 10)

;; Succ edge  10 [100.0%] 

(barrier 44 43 45)

;; Start of basic block ( 6) -> 8
;; Pred edge  6 [95.5%] 
(code_label 45 44 46 8 96 "" [1 uses])

(note 46 45 47 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 47 46 48 8 sequenceUtil.cpp:32 (parallel [
            (set (reg/v:SI 61 [ i ])
                (plus:SI (reg/v:SI 61 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(debug_insn 48 47 49 8 sequenceUtil.cpp:32 (var_location:SI i (reg/v:SI 61 [ i ])) -1 (nil))

(debug_insn 49 48 50 8 (var_location:SI i (reg/v:SI 61 [ i ])) -1 (nil))

(insn 50 49 51 8 sequenceUtil.cpp:32 (set (reg:DI 63 [ D.31789 ])
        (zero_extend:DI (reg/v:SI 61 [ i ]))) -1 (nil))

(debug_insn 51 50 52 8 sequenceUtil.cpp:32 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 52 51 53 8 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(debug_insn 53 52 54 8 (var_location:DI this (reg/v/f:DI 65 [ kmer ])) -1 (nil))

(insn 54 53 56 8 /usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../include/c++/4.4.6/bits/basic_string.h:286 (parallel [
            (set (reg/f:DI 60 [ D.39717 ])
                (plus:DI (reg/f:DI 59 [ prephitmp.887 ])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil))

(insn 56 54 57 8 sequenceUtil.cpp:32 (set (reg:CC 17 flags)
        (compare:CC (reg:DI 63 [ D.31789 ])
            (mem/s:DI (reg/f:DI 60 [ D.39717 ]) [14 <variable>.D.11486._M_length+0 S8 A64]))) -1 (nil))

(jump_insn 57 56 58 8 sequenceUtil.cpp:32 (set (pc)
        (if_then_else (ltu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 55)
            (pc))) -1 (expr_list:REG_BR_PROB (const_int 9550 [0x254e])
        (nil)))
;; End of basic block 8 -> ( 4 9)

;; Succ edge  4 [95.5%]  (dfs_back)
;; Succ edge  9 [4.5%]  (fallthru)

;; Start of basic block ( 2 8) -> 9
;; Pred edge  2 [4.5%] 
;; Pred edge  8 [4.5%]  (fallthru)
(code_label 58 57 59 9 94 "" [1 uses])

(note 59 58 60 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 60 59 61 9 sequenceUtil.cpp:32 (set (reg:QI 62 [ D.31798 ])
        (const_int 0 [0x0])) -1 (nil))
;; End of basic block 9 -> ( 10)

;; Succ edge  10 [100.0%]  (fallthru)

;; Start of basic block ( 9 7) -> 10
;; Pred edge  9 [100.0%]  (fallthru)
;; Pred edge  7 [100.0%] 
(code_label 61 60 62 10 97 "" [1 uses])

(note 62 61 63 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(insn 63 62 67 10 sequenceUtil.cpp:32 (set (reg:QI 64 [ <result> ])
        (reg:QI 62 [ D.31798 ])) -1 (nil))

(insn 67 63 73 10 sequenceUtil.cpp:49 (set (reg/i:QI 0 ax)
        (reg:QI 64 [ <result> ])) -1 (nil))

(insn 73 67 0 10 sequenceUtil.cpp:49 (use (reg/i:QI 0 ax)) -1 (nil))
;; End of basic block 10 -> ( 1)

;; Succ edge  EXIT [100.0%]  (fallthru)
