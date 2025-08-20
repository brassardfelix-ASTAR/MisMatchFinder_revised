; ModuleID = 'probe1.47fb7fe357dd304b-cgu.0'
source_filename = "probe1.47fb7fe357dd304b-cgu.0"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128"
target triple = "x86_64-unknown-linux-gnu"

%"core::fmt::rt::Argument<'_>" = type { %"core::fmt::rt::ArgumentType<'_>" }
%"core::fmt::rt::ArgumentType<'_>" = type { ptr, [1 x i64] }

@anon.82179aa89ed86be757434e5965125d80.0 = private unnamed_addr constant <{ [8 x i8], [8 x i8] }> <{ [8 x i8] zeroinitializer, [8 x i8] undef }>, align 8
@alloc_fad0cd83b7d1858a846a172eb260e593 = private unnamed_addr constant [42 x i8] c"is_aligned_to: align is not a power-of-two", align 1
@alloc_e92e94d0ff530782b571cfd99ec66aef = private unnamed_addr constant <{ ptr, [8 x i8] }> <{ ptr @alloc_fad0cd83b7d1858a846a172eb260e593, [8 x i8] c"*\00\00\00\00\00\00\00" }>, align 8
@alloc_a91c3d56db4a4aea2768c8f12b410bbd = private unnamed_addr constant [82 x i8] c"/rustc/8e62bfd311791bfd9dca886abdfbab07ec54d8b4/library/core/src/ptr/const_ptr.rs\00", align 1
@alloc_59d15bf41af9dd082a1bb68a07ee1f5c = private unnamed_addr constant <{ ptr, [16 x i8] }> <{ ptr @alloc_a91c3d56db4a4aea2768c8f12b410bbd, [16 x i8] c"Q\00\00\00\00\00\00\00\92\05\00\00\0D\00\00\00" }>, align 8
@alloc_bd3468a7b96187f70c1ce98a3e7a63bf = private unnamed_addr constant [283 x i8] c"unsafe precondition(s) violated: ptr::copy_nonoverlapping requires that both pointer arguments are aligned and non-null and the specified memory ranges do not overlap\0A\0AThis indicates a bug in the program. This Undefined Behavior check is optional, and cannot be relied on for safety.", align 1
@alloc_64e308ef4babfeb8b6220184de794a17 = private unnamed_addr constant [221 x i8] c"unsafe precondition(s) violated: hint::assert_unchecked must never be called when the condition is false\0A\0AThis indicates a bug in the program. This Undefined Behavior check is optional, and cannot be relied on for safety.", align 1
@alloc_763310d78c99c2c1ad3f8a9821e942f3 = private unnamed_addr constant [61 x i8] c"is_nonoverlapping: `size_of::<T>() * count` overflows a usize", align 1
@alloc_57d70e9d94c65ecfc15225d29a5ed72b = private unnamed_addr constant [198 x i8] c"unsafe precondition(s) violated: Vec::set_len requires that new_len <= capacity()\0A\0AThis indicates a bug in the program. This Undefined Behavior check is optional, and cannot be relied on for safety.", align 1
@alloc_b6c5d1263b121da531320d7d71ce14dc = private unnamed_addr constant [75 x i8] c"/rustc/8e62bfd311791bfd9dca886abdfbab07ec54d8b4/library/alloc/src/slice.rs\00", align 1
@alloc_ac65f47950d5bfcf522958a117ed0fd0 = private unnamed_addr constant <{ ptr, [16 x i8] }> <{ ptr @alloc_b6c5d1263b121da531320d7d71ce14dc, [16 x i8] c"J\00\00\00\00\00\00\00\BD\01\00\00\1D\00\00\00" }>, align 8
@alloc_73d79880c03600f0590ae12b69e696d2 = private unnamed_addr constant [76 x i8] c"/rustc/8e62bfd311791bfd9dca886abdfbab07ec54d8b4/library/core/src/ptr/mod.rs\00", align 1
@alloc_b72098898095f52adc73e53a589e397d = private unnamed_addr constant <{ ptr, [16 x i8] }> <{ ptr @alloc_73d79880c03600f0590ae12b69e696d2, [16 x i8] c"K\00\00\00\00\00\00\00\0F\02\00\00\05\00\00\00" }>, align 8
@alloc_6bd2faa7b3a88110100b00c9c177bd0f = private unnamed_addr constant [77 x i8] c"/rustc/8e62bfd311791bfd9dca886abdfbab07ec54d8b4/library/alloc/src/vec/mod.rs\00", align 1
@alloc_2777e455911d18e2781b88be07504468 = private unnamed_addr constant <{ ptr, [16 x i8] }> <{ ptr @alloc_6bd2faa7b3a88110100b00c9c177bd0f, [16 x i8] c"L\00\00\00\00\00\00\00\DD\07\00\00\09\00\00\00" }>, align 8
@alloc_53973d2fe29b4adba8bb7390b5678745 = private unnamed_addr constant [8 x i8] zeroinitializer, align 8
@alloc_eab5d04767146d7d9b93b60d28ef530a = private unnamed_addr constant <{ ptr, [8 x i8] }> <{ ptr inttoptr (i64 1 to ptr), [8 x i8] zeroinitializer }>, align 8

; core::intrinsics::cold_path
; Function Attrs: cold nounwind nonlazybind uwtable
define internal void @_ZN4core10intrinsics9cold_path17ha3eab03867a0cc13E() unnamed_addr #0 {
start:
  ret void
}

; core::fmt::rt::<impl core::fmt::Arguments>::new_v1
; Function Attrs: inlinehint nonlazybind uwtable
define void @"_ZN4core3fmt2rt38_$LT$impl$u20$core..fmt..Arguments$GT$6new_v117hcf32058b4dd71b39E"(ptr sret([48 x i8]) align 8 %_0, ptr align 8 %pieces, ptr align 8 %args) unnamed_addr #1 {
start:
  store ptr %pieces, ptr %_0, align 8
  %0 = getelementptr inbounds i8, ptr %_0, i64 8
  store i64 1, ptr %0, align 8
  %1 = load ptr, ptr @anon.82179aa89ed86be757434e5965125d80.0, align 8
  %2 = load i64, ptr getelementptr inbounds (i8, ptr @anon.82179aa89ed86be757434e5965125d80.0, i64 8), align 8
  %3 = getelementptr inbounds i8, ptr %_0, i64 32
  store ptr %1, ptr %3, align 8
  %4 = getelementptr inbounds i8, ptr %3, i64 8
  store i64 %2, ptr %4, align 8
  %5 = getelementptr inbounds i8, ptr %_0, i64 16
  store ptr %args, ptr %5, align 8
  %6 = getelementptr inbounds i8, ptr %5, i64 8
  store i64 1, ptr %6, align 8
  ret void
}

; core::fmt::rt::Argument::new_lower_exp
; Function Attrs: inlinehint nonlazybind uwtable
define void @_ZN4core3fmt2rt8Argument13new_lower_exp17hdb28f0a75dad489fE(ptr sret([16 x i8]) align 8 %_0, ptr align 8 %x) unnamed_addr #1 {
start:
  %_2 = alloca [16 x i8], align 8
  store ptr %x, ptr %_2, align 8
  %0 = getelementptr inbounds i8, ptr %_2, i64 8
  store ptr @"_ZN4core3fmt3num3imp55_$LT$impl$u20$core..fmt..LowerExp$u20$for$u20$isize$GT$3fmt17h47b185332d7ed1c2E", ptr %0, align 8
  call void @llvm.memcpy.p0.p0.i64(ptr align 8 %_0, ptr align 8 %_2, i64 16, i1 false)
  ret void
}

; core::ops::function::FnOnce::call_once
; Function Attrs: inlinehint nonlazybind uwtable
define internal void @_ZN4core3ops8function6FnOnce9call_once17h3f0767c35fad6ab0E(ptr sret([24 x i8]) align 8 %_0, ptr align 1 %0, i64 %1) unnamed_addr #1 {
start:
  %_2 = alloca [16 x i8], align 8
  store ptr %0, ptr %_2, align 8
  %2 = getelementptr inbounds i8, ptr %_2, i64 8
  store i64 %1, ptr %2, align 8
  %3 = load ptr, ptr %_2, align 8
  %4 = getelementptr inbounds i8, ptr %_2, i64 8
  %5 = load i64, ptr %4, align 8
; call alloc::str::<impl alloc::borrow::ToOwned for str>::to_owned
  call void @"_ZN5alloc3str56_$LT$impl$u20$alloc..borrow..ToOwned$u20$for$u20$str$GT$8to_owned17ha410b0efbfe4e8fcE"(ptr sret([24 x i8]) align 8 %_0, ptr align 1 %3, i64 %5)
  ret void
}

; core::ptr::copy_nonoverlapping::precondition_check
; Function Attrs: inlinehint nounwind nonlazybind uwtable
define internal void @_ZN4core3ptr19copy_nonoverlapping18precondition_check17hd1660667b8df50cdE(ptr %src, ptr %dst, i64 %size, i64 %align, i64 %count, ptr align 8 %0) unnamed_addr #2 personality ptr @rust_eh_personality {
start:
  %1 = alloca [4 x i8], align 4
  %_29 = alloca [48 x i8], align 8
  %_24 = alloca [4 x i8], align 4
  %_23 = alloca [8 x i8], align 8
  %_22 = alloca [8 x i8], align 8
  %_21 = alloca [8 x i8], align 8
  %_20 = alloca [48 x i8], align 8
  %_17 = alloca [16 x i8], align 8
  %_15 = alloca [48 x i8], align 8
  %is_zst = alloca [1 x i8], align 1
  %align1 = alloca [8 x i8], align 8
  %zero_size = alloca [1 x i8], align 1
  %2 = icmp eq i64 %count, 0
  br i1 %2, label %bb1, label %bb2

bb1:                                              ; preds = %start
  store i8 1, ptr %zero_size, align 1
  store i64 %align, ptr %align1, align 8
  %3 = load i8, ptr %zero_size, align 1
  %4 = trunc nuw i8 %3 to i1
  %5 = zext i1 %4 to i8
  store i8 %5, ptr %is_zst, align 1
  %6 = call i64 @llvm.ctpop.i64(i64 %align)
  %7 = trunc i64 %6 to i32
  store i32 %7, ptr %_24, align 4
  %8 = load i32, ptr %_24, align 4
  %9 = icmp eq i32 %8, 1
  br i1 %9, label %bb26, label %bb15

bb2:                                              ; preds = %start
  %10 = icmp eq i64 %size, 0
  %11 = zext i1 %10 to i8
  store i8 %11, ptr %zero_size, align 1
  store i64 %align, ptr %align1, align 8
  %12 = load i8, ptr %zero_size, align 1
  %13 = trunc nuw i8 %12 to i1
  %14 = zext i1 %13 to i8
  store i8 %14, ptr %is_zst, align 1
  %15 = call i64 @llvm.ctpop.i64(i64 %align)
  %16 = trunc i64 %15 to i32
  store i32 %16, ptr %_24, align 4
  %17 = load i32, ptr %_24, align 4
  %18 = icmp eq i32 %17, 1
  br i1 %18, label %bb14, label %bb15

bb26:                                             ; preds = %bb1
  store ptr %src, ptr %_22, align 8
  %19 = sub i64 %align, 1
  store i64 %19, ptr %_23, align 8
  %20 = load i64, ptr %_22, align 8
  %21 = load i64, ptr %_23, align 8
  %22 = and i64 %20, %21
  store i64 %22, ptr %_21, align 8
  %23 = load i64, ptr %_21, align 8
  %24 = icmp eq i64 %23, 0
  br i1 %24, label %bb27, label %bb11

bb15:                                             ; preds = %bb2, %bb1
  store ptr @alloc_e92e94d0ff530782b571cfd99ec66aef, ptr %_20, align 8
  %25 = getelementptr inbounds i8, ptr %_20, i64 8
  store i64 1, ptr %25, align 8
  %26 = load ptr, ptr @anon.82179aa89ed86be757434e5965125d80.0, align 8
  %27 = load i64, ptr getelementptr inbounds (i8, ptr @anon.82179aa89ed86be757434e5965125d80.0, i64 8), align 8
  %28 = getelementptr inbounds i8, ptr %_20, i64 32
  store ptr %26, ptr %28, align 8
  %29 = getelementptr inbounds i8, ptr %28, i64 8
  store i64 %27, ptr %29, align 8
  %30 = getelementptr inbounds i8, ptr %_20, i64 16
  store ptr inttoptr (i64 8 to ptr), ptr %30, align 8
  %31 = getelementptr inbounds i8, ptr %30, i64 8
  store i64 0, ptr %31, align 8
; invoke core::panicking::panic_fmt
  invoke void @_ZN4core9panicking9panic_fmt17h0023f2c1f721ac1dE(ptr align 8 %_20, ptr align 8 @alloc_59d15bf41af9dd082a1bb68a07ee1f5c) #11
          to label %unreachable unwind label %terminate

bb27:                                             ; preds = %bb26
  br label %bb12

bb11:                                             ; preds = %bb14, %bb26
  br label %bb6

bb12:                                             ; preds = %bb10, %bb27
  br label %bb3

bb14:                                             ; preds = %bb2
  store ptr %src, ptr %_22, align 8
  %32 = sub i64 %align, 1
  store i64 %32, ptr %_23, align 8
  %33 = load i64, ptr %_22, align 8
  %34 = load i64, ptr %_23, align 8
  %35 = and i64 %33, %34
  store i64 %35, ptr %_21, align 8
  %36 = load i64, ptr %_21, align 8
  %37 = icmp eq i64 %36, 0
  br i1 %37, label %bb10, label %bb11

bb10:                                             ; preds = %bb14
  %38 = load i8, ptr %is_zst, align 1
  %39 = trunc nuw i8 %38 to i1
  br i1 %39, label %bb12, label %bb13

bb13:                                             ; preds = %bb10
  %40 = load i64, ptr %_22, align 8
  %_18 = icmp eq i64 %40, 0
  %_8 = xor i1 %_18, true
  br i1 %_8, label %bb3, label %bb6

bb6:                                              ; preds = %bb11, %bb13
  br label %bb7

bb3:                                              ; preds = %bb12, %bb13
  %41 = load i8, ptr %zero_size, align 1
  %is_zst2 = trunc nuw i8 %41 to i1
  %42 = call i64 @llvm.ctpop.i64(i64 %align)
  %43 = trunc i64 %42 to i32
  store i32 %43, ptr %1, align 4
  %_32 = load i32, ptr %1, align 4
  %44 = icmp eq i32 %_32, 1
  br i1 %44, label %bb21, label %bb22

bb21:                                             ; preds = %bb3
  %_31 = ptrtoint ptr %dst to i64
  %45 = load i64, ptr %_23, align 8
  %_30 = and i64 %_31, %45
  %46 = icmp eq i64 %_30, 0
  br i1 %46, label %bb17, label %bb18

bb22:                                             ; preds = %bb3
  store ptr @alloc_e92e94d0ff530782b571cfd99ec66aef, ptr %_29, align 8
  %47 = getelementptr inbounds i8, ptr %_29, i64 8
  store i64 1, ptr %47, align 8
  %48 = load ptr, ptr @anon.82179aa89ed86be757434e5965125d80.0, align 8
  %49 = load i64, ptr getelementptr inbounds (i8, ptr @anon.82179aa89ed86be757434e5965125d80.0, i64 8), align 8
  %50 = getelementptr inbounds i8, ptr %_29, i64 32
  store ptr %48, ptr %50, align 8
  %51 = getelementptr inbounds i8, ptr %50, i64 8
  store i64 %49, ptr %51, align 8
  %52 = getelementptr inbounds i8, ptr %_29, i64 16
  store ptr inttoptr (i64 8 to ptr), ptr %52, align 8
  %53 = getelementptr inbounds i8, ptr %52, i64 8
  store i64 0, ptr %53, align 8
; invoke core::panicking::panic_fmt
  invoke void @_ZN4core9panicking9panic_fmt17h0023f2c1f721ac1dE(ptr align 8 %_29, ptr align 8 @alloc_59d15bf41af9dd082a1bb68a07ee1f5c) #11
          to label %unreachable unwind label %terminate

bb17:                                             ; preds = %bb21
  br i1 %is_zst2, label %bb19, label %bb20

bb18:                                             ; preds = %bb21
  br label %bb5

bb20:                                             ; preds = %bb17
  %_27 = icmp eq i64 %_31, 0
  %_11 = xor i1 %_27, true
  br i1 %_11, label %bb4, label %bb5

bb19:                                             ; preds = %bb17
  br label %bb4

bb5:                                              ; preds = %bb18, %bb20
  br label %bb7

bb4:                                              ; preds = %bb19, %bb20
; invoke core::ub_checks::maybe_is_nonoverlapping::runtime
  %_6 = invoke zeroext i1 @_ZN4core9ub_checks23maybe_is_nonoverlapping7runtime17h625bc961682578f4E(ptr %src, ptr %dst, i64 %size, i64 %count)
          to label %bb24 unwind label %terminate

terminate:                                        ; preds = %bb15, %bb22, %bb4
  %54 = landingpad { ptr, i32 }
          filter [0 x ptr] zeroinitializer
; call core::panicking::panic_cannot_unwind
  call void @_ZN4core9panicking19panic_cannot_unwind17h48697f70b892f57bE() #12
  unreachable

bb24:                                             ; preds = %bb4
  br i1 %_6, label %bb9, label %bb8

bb8:                                              ; preds = %bb7, %bb24
  %55 = getelementptr inbounds nuw { ptr, i64 }, ptr %_17, i64 0
  store ptr @alloc_bd3468a7b96187f70c1ce98a3e7a63bf, ptr %55, align 8
  %56 = getelementptr inbounds i8, ptr %55, i64 8
  store i64 283, ptr %56, align 8
  store ptr %_17, ptr %_15, align 8
  %57 = getelementptr inbounds i8, ptr %_15, i64 8
  store i64 1, ptr %57, align 8
  %58 = load ptr, ptr @anon.82179aa89ed86be757434e5965125d80.0, align 8
  %59 = load i64, ptr getelementptr inbounds (i8, ptr @anon.82179aa89ed86be757434e5965125d80.0, i64 8), align 8
  %60 = getelementptr inbounds i8, ptr %_15, i64 32
  store ptr %58, ptr %60, align 8
  %61 = getelementptr inbounds i8, ptr %60, i64 8
  store i64 %59, ptr %61, align 8
  %62 = getelementptr inbounds i8, ptr %_15, i64 16
  store ptr inttoptr (i64 8 to ptr), ptr %62, align 8
  %63 = getelementptr inbounds i8, ptr %62, i64 8
  store i64 0, ptr %63, align 8
; call core::panicking::panic_nounwind_fmt
  call void @_ZN4core9panicking18panic_nounwind_fmt17h219ec714db0583beE(ptr align 8 %_15, i1 zeroext false, ptr align 8 %0) #13
  unreachable

bb9:                                              ; preds = %bb24
  ret void

bb7:                                              ; preds = %bb6, %bb5
  br label %bb8

unreachable:                                      ; preds = %bb15, %bb22
  unreachable
}

; core::ptr::drop_in_place<alloc::string::String>
; Function Attrs: nonlazybind uwtable
define void @"_ZN4core3ptr42drop_in_place$LT$alloc..string..String$GT$17hacbd5aff21b03340E"(ptr align 8 %_1) unnamed_addr #3 {
start:
; call core::ptr::drop_in_place<alloc::vec::Vec<u8>>
  call void @"_ZN4core3ptr46drop_in_place$LT$alloc..vec..Vec$LT$u8$GT$$GT$17hd9985c0cac6ae0e5E"(ptr align 8 %_1)
  ret void
}

; core::ptr::drop_in_place<alloc::vec::Vec<u8>>
; Function Attrs: nonlazybind uwtable
define void @"_ZN4core3ptr46drop_in_place$LT$alloc..vec..Vec$LT$u8$GT$$GT$17hd9985c0cac6ae0e5E"(ptr align 8 %_1) unnamed_addr #3 personality ptr @rust_eh_personality {
start:
  %0 = alloca [16 x i8], align 8
; invoke <alloc::vec::Vec<T,A> as core::ops::drop::Drop>::drop
  invoke void @"_ZN70_$LT$alloc..vec..Vec$LT$T$C$A$GT$$u20$as$u20$core..ops..drop..Drop$GT$4drop17h991b5fea7f8f6006E"(ptr align 8 %_1)
          to label %bb4 unwind label %cleanup

bb3:                                              ; preds = %cleanup
; invoke core::ptr::drop_in_place<alloc::raw_vec::RawVec<u8>>
  invoke void @"_ZN4core3ptr53drop_in_place$LT$alloc..raw_vec..RawVec$LT$u8$GT$$GT$17ha1ffb92fde325907E"(ptr align 8 %_1) #14
          to label %bb1 unwind label %terminate

cleanup:                                          ; preds = %start
  %1 = landingpad { ptr, i32 }
          cleanup
  %2 = extractvalue { ptr, i32 } %1, 0
  %3 = extractvalue { ptr, i32 } %1, 1
  store ptr %2, ptr %0, align 8
  %4 = getelementptr inbounds i8, ptr %0, i64 8
  store i32 %3, ptr %4, align 8
  br label %bb3

bb4:                                              ; preds = %start
; call core::ptr::drop_in_place<alloc::raw_vec::RawVec<u8>>
  call void @"_ZN4core3ptr53drop_in_place$LT$alloc..raw_vec..RawVec$LT$u8$GT$$GT$17ha1ffb92fde325907E"(ptr align 8 %_1)
  ret void

terminate:                                        ; preds = %bb3
  %5 = landingpad { ptr, i32 }
          filter [0 x ptr] zeroinitializer
; call core::panicking::panic_in_cleanup
  call void @_ZN4core9panicking16panic_in_cleanup17h7a69428c8f5fb9a4E() #12
  unreachable

bb1:                                              ; preds = %bb3
  %6 = load ptr, ptr %0, align 8
  %7 = getelementptr inbounds i8, ptr %0, i64 8
  %8 = load i32, ptr %7, align 8
  %9 = insertvalue { ptr, i32 } poison, ptr %6, 0
  %10 = insertvalue { ptr, i32 } %9, i32 %8, 1
  resume { ptr, i32 } %10
}

; core::ptr::drop_in_place<alloc::raw_vec::RawVec<u8>>
; Function Attrs: nonlazybind uwtable
define void @"_ZN4core3ptr53drop_in_place$LT$alloc..raw_vec..RawVec$LT$u8$GT$$GT$17ha1ffb92fde325907E"(ptr align 8 %_1) unnamed_addr #3 {
start:
; call <alloc::raw_vec::RawVec<T,A> as core::ops::drop::Drop>::drop
  call void @"_ZN77_$LT$alloc..raw_vec..RawVec$LT$T$C$A$GT$$u20$as$u20$core..ops..drop..Drop$GT$4drop17h902845a8e4723e86E"(ptr align 8 %_1)
  ret void
}

; core::hint::assert_unchecked::precondition_check
; Function Attrs: inlinehint nounwind nonlazybind uwtable
define internal void @_ZN4core4hint16assert_unchecked18precondition_check17h979bd9711919bbd4E(i1 zeroext %cond, ptr align 8 %0) unnamed_addr #2 {
start:
  %_5 = alloca [16 x i8], align 8
  %_3 = alloca [48 x i8], align 8
  br i1 %cond, label %bb2, label %bb1

bb1:                                              ; preds = %start
  %1 = getelementptr inbounds nuw { ptr, i64 }, ptr %_5, i64 0
  store ptr @alloc_64e308ef4babfeb8b6220184de794a17, ptr %1, align 8
  %2 = getelementptr inbounds i8, ptr %1, i64 8
  store i64 221, ptr %2, align 8
  store ptr %_5, ptr %_3, align 8
  %3 = getelementptr inbounds i8, ptr %_3, i64 8
  store i64 1, ptr %3, align 8
  %4 = load ptr, ptr @anon.82179aa89ed86be757434e5965125d80.0, align 8
  %5 = load i64, ptr getelementptr inbounds (i8, ptr @anon.82179aa89ed86be757434e5965125d80.0, i64 8), align 8
  %6 = getelementptr inbounds i8, ptr %_3, i64 32
  store ptr %4, ptr %6, align 8
  %7 = getelementptr inbounds i8, ptr %6, i64 8
  store i64 %5, ptr %7, align 8
  %8 = getelementptr inbounds i8, ptr %_3, i64 16
  store ptr inttoptr (i64 8 to ptr), ptr %8, align 8
  %9 = getelementptr inbounds i8, ptr %8, i64 8
  store i64 0, ptr %9, align 8
; call core::panicking::panic_nounwind_fmt
  call void @_ZN4core9panicking18panic_nounwind_fmt17h219ec714db0583beE(ptr align 8 %_3, i1 zeroext false, ptr align 8 %0) #13
  unreachable

bb2:                                              ; preds = %start
  ret void
}

; core::option::Option<T>::map_or_else
; Function Attrs: inlinehint nonlazybind uwtable
define void @"_ZN4core6option15Option$LT$T$GT$11map_or_else17h37ce12d1eeecdc62E"(ptr sret([24 x i8]) align 8 %_0, ptr align 1 %0, i64 %1, ptr align 8 %default) unnamed_addr #1 personality ptr @rust_eh_personality {
start:
  %2 = alloca [16 x i8], align 8
  %_10 = alloca [1 x i8], align 1
  %_9 = alloca [1 x i8], align 1
  %self = alloca [16 x i8], align 8
  store ptr %0, ptr %self, align 8
  %3 = getelementptr inbounds i8, ptr %self, i64 8
  store i64 %1, ptr %3, align 8
  store i8 1, ptr %_10, align 1
  store i8 1, ptr %_9, align 1
  %4 = load ptr, ptr %self, align 8
  %5 = getelementptr inbounds i8, ptr %self, i64 8
  %6 = load i64, ptr %5, align 8
  %7 = ptrtoint ptr %4 to i64
  %8 = icmp eq i64 %7, 0
  %_4 = select i1 %8, i64 0, i64 1
  %9 = trunc nuw i64 %_4 to i1
  br i1 %9, label %bb3, label %bb2

bb3:                                              ; preds = %start
  %t.0 = load ptr, ptr %self, align 8
  %10 = getelementptr inbounds i8, ptr %self, i64 8
  %t.1 = load i64, ptr %10, align 8
  store i8 0, ptr %_9, align 1
; invoke core::ops::function::FnOnce::call_once
  invoke void @_ZN4core3ops8function6FnOnce9call_once17h3f0767c35fad6ab0E(ptr sret([24 x i8]) align 8 %_0, ptr align 1 %t.0, i64 %t.1)
          to label %bb4 unwind label %cleanup

bb2:                                              ; preds = %start
  store i8 0, ptr %_10, align 1
; invoke alloc::fmt::format::{{closure}}
  invoke void @"_ZN5alloc3fmt6format28_$u7b$$u7b$closure$u7d$$u7d$17h38a57ff6c8acb163E"(ptr sret([24 x i8]) align 8 %_0, ptr align 8 %default)
          to label %bb5 unwind label %cleanup

bb11:                                             ; preds = %cleanup
  %11 = load i8, ptr %_9, align 1
  %12 = trunc nuw i8 %11 to i1
  br i1 %12, label %bb10, label %bb7

cleanup:                                          ; preds = %bb3, %bb2
  %13 = landingpad { ptr, i32 }
          cleanup
  %14 = extractvalue { ptr, i32 } %13, 0
  %15 = extractvalue { ptr, i32 } %13, 1
  store ptr %14, ptr %2, align 8
  %16 = getelementptr inbounds i8, ptr %2, i64 8
  store i32 %15, ptr %16, align 8
  br label %bb11

bb5:                                              ; preds = %bb2
  br label %bb6

bb6:                                              ; preds = %bb9, %bb4, %bb5
  ret void

bb4:                                              ; preds = %bb3
  %17 = load i8, ptr %_10, align 1
  %18 = trunc nuw i8 %17 to i1
  br i1 %18, label %bb9, label %bb6

bb9:                                              ; preds = %bb4
  br label %bb6

bb7:                                              ; preds = %bb10, %bb11
  %19 = load i8, ptr %_10, align 1
  %20 = trunc nuw i8 %19 to i1
  br i1 %20, label %bb12, label %bb8

bb10:                                             ; preds = %bb11
  br label %bb7

bb8:                                              ; preds = %bb12, %bb7
  %21 = load ptr, ptr %2, align 8
  %22 = getelementptr inbounds i8, ptr %2, i64 8
  %23 = load i32, ptr %22, align 8
  %24 = insertvalue { ptr, i32 } poison, ptr %21, 0
  %25 = insertvalue { ptr, i32 } %24, i32 %23, 1
  resume { ptr, i32 } %25

bb12:                                             ; preds = %bb7
  br label %bb8

bb1:                                              ; No predecessors!
  unreachable
}

; core::ub_checks::maybe_is_nonoverlapping::runtime
; Function Attrs: inlinehint nonlazybind uwtable
define internal zeroext i1 @_ZN4core9ub_checks23maybe_is_nonoverlapping7runtime17h625bc961682578f4E(ptr %src, ptr %dst, i64 %size, i64 %count) unnamed_addr #1 {
start:
  %diff = alloca [8 x i8], align 8
  %_9 = alloca [16 x i8], align 8
  %src_usize = ptrtoint ptr %src to i64
  %dst_usize = ptrtoint ptr %dst to i64
  %0 = call { i64, i1 } @llvm.umul.with.overflow.i64(i64 %size, i64 %count)
  %_14.0 = extractvalue { i64, i1 } %0, 0
  %_14.1 = extractvalue { i64, i1 } %0, 1
  br i1 %_14.1, label %bb1, label %bb3

bb3:                                              ; preds = %start
  %1 = getelementptr inbounds i8, ptr %_9, i64 8
  store i64 %_14.0, ptr %1, align 8
  store i64 1, ptr %_9, align 8
  %2 = getelementptr inbounds i8, ptr %_9, i64 8
  %size1 = load i64, ptr %2, align 8
  %_22 = icmp ult i64 %src_usize, %dst_usize
  br i1 %_22, label %bb4, label %bb5

bb1:                                              ; preds = %start
; call core::panicking::panic_nounwind
  call void @_ZN4core9panicking14panic_nounwind17h47fee6258c2f7571E(ptr align 1 @alloc_763310d78c99c2c1ad3f8a9821e942f3, i64 61) #13
  unreachable

bb5:                                              ; preds = %bb3
  %3 = sub i64 %src_usize, %dst_usize
  store i64 %3, ptr %diff, align 8
  br label %bb6

bb4:                                              ; preds = %bb3
  %4 = sub i64 %dst_usize, %src_usize
  store i64 %4, ptr %diff, align 8
  br label %bb6

bb6:                                              ; preds = %bb4, %bb5
  %_11 = load i64, ptr %diff, align 8
  %_0 = icmp uge i64 %_11, %size1
  ret i1 %_0
}

; alloc::fmt::format
; Function Attrs: inlinehint nonlazybind uwtable
define internal void @_ZN5alloc3fmt6format17h6ed158a03a5de8a7E(ptr sret([24 x i8]) align 8 %_0, ptr align 8 %args) unnamed_addr #1 {
start:
  %_2 = alloca [16 x i8], align 8
  %_6.0 = load ptr, ptr %args, align 8
  %0 = getelementptr inbounds i8, ptr %args, i64 8
  %_6.1 = load i64, ptr %0, align 8
  %1 = getelementptr inbounds i8, ptr %args, i64 16
  %_7.0 = load ptr, ptr %1, align 8
  %2 = getelementptr inbounds i8, ptr %1, i64 8
  %_7.1 = load i64, ptr %2, align 8
  %3 = icmp eq i64 %_6.1, 0
  br i1 %3, label %bb4, label %bb5

bb4:                                              ; preds = %start
  %4 = icmp eq i64 %_7.1, 0
  br i1 %4, label %bb8, label %bb3

bb5:                                              ; preds = %start
  %5 = icmp eq i64 %_6.1, 1
  br i1 %5, label %bb6, label %bb3

bb8:                                              ; preds = %bb4
  store ptr inttoptr (i64 1 to ptr), ptr %_2, align 8
  %6 = getelementptr inbounds i8, ptr %_2, i64 8
  store i64 0, ptr %6, align 8
  br label %bb2

bb3:                                              ; preds = %bb6, %bb5, %bb4
  %7 = load ptr, ptr @anon.82179aa89ed86be757434e5965125d80.0, align 8
  %8 = load i64, ptr getelementptr inbounds (i8, ptr @anon.82179aa89ed86be757434e5965125d80.0, i64 8), align 8
  store ptr %7, ptr %_2, align 8
  %9 = getelementptr inbounds i8, ptr %_2, i64 8
  store i64 %8, ptr %9, align 8
  br label %bb2

bb2:                                              ; preds = %bb3, %bb7, %bb8
  %10 = load ptr, ptr %_2, align 8
  %11 = getelementptr inbounds i8, ptr %_2, i64 8
  %12 = load i64, ptr %11, align 8
; call core::option::Option<T>::map_or_else
  call void @"_ZN4core6option15Option$LT$T$GT$11map_or_else17h37ce12d1eeecdc62E"(ptr sret([24 x i8]) align 8 %_0, ptr align 1 %10, i64 %12, ptr align 8 %args)
  ret void

bb6:                                              ; preds = %bb5
  %13 = icmp eq i64 %_7.1, 0
  br i1 %13, label %bb7, label %bb3

bb7:                                              ; preds = %bb6
  %s = getelementptr inbounds nuw { ptr, i64 }, ptr %_6.0, i64 0
  %14 = getelementptr inbounds nuw { ptr, i64 }, ptr %_6.0, i64 0
  %_12.0 = load ptr, ptr %14, align 8
  %15 = getelementptr inbounds i8, ptr %14, i64 8
  %_12.1 = load i64, ptr %15, align 8
  store ptr %_12.0, ptr %_2, align 8
  %16 = getelementptr inbounds i8, ptr %_2, i64 8
  store i64 %_12.1, ptr %16, align 8
  br label %bb2
}

; alloc::fmt::format::{{closure}}
; Function Attrs: inlinehint nonlazybind uwtable
define void @"_ZN5alloc3fmt6format28_$u7b$$u7b$closure$u7d$$u7d$17h38a57ff6c8acb163E"(ptr sret([24 x i8]) align 8 %_0, ptr align 8 %_1) unnamed_addr #1 {
start:
  %_2 = alloca [48 x i8], align 8
  call void @llvm.memcpy.p0.p0.i64(ptr align 8 %_2, ptr align 8 %_1, i64 48, i1 false)
; call alloc::fmt::format::format_inner
  call void @_ZN5alloc3fmt6format12format_inner17h52d1570151877af0E(ptr sret([24 x i8]) align 8 %_0, ptr align 8 %_2)
  ret void
}

; alloc::str::<impl alloc::borrow::ToOwned for str>::to_owned
; Function Attrs: inlinehint nonlazybind uwtable
define internal void @"_ZN5alloc3str56_$LT$impl$u20$alloc..borrow..ToOwned$u20$for$u20$str$GT$8to_owned17ha410b0efbfe4e8fcE"(ptr sret([24 x i8]) align 8 %_0, ptr align 1 %self.0, i64 %self.1) unnamed_addr #1 {
start:
  %bytes = alloca [24 x i8], align 8
; call <T as alloc::slice::<impl [T]>::to_vec_in::ConvertVec>::to_vec
  call void @"_ZN87_$LT$T$u20$as$u20$alloc..slice..$LT$impl$u20$$u5b$T$u5d$$GT$..to_vec_in..ConvertVec$GT$6to_vec17h6b2167f36aa25ff3E"(ptr sret([24 x i8]) align 8 %bytes, ptr align 1 %self.0, i64 %self.1)
  call void @llvm.memcpy.p0.p0.i64(ptr align 8 %_0, ptr align 8 %bytes, i64 24, i1 false)
  ret void
}

; alloc::vec::Vec<T,A>::set_len::precondition_check
; Function Attrs: inlinehint nounwind nonlazybind uwtable
define internal void @"_ZN5alloc3vec16Vec$LT$T$C$A$GT$7set_len18precondition_check17h36c85061f27ed97bE"(i64 %new_len, i64 %capacity, ptr align 8 %0) unnamed_addr #2 {
start:
  %_7 = alloca [16 x i8], align 8
  %_5 = alloca [48 x i8], align 8
  %_3 = icmp ule i64 %new_len, %capacity
  br i1 %_3, label %bb1, label %bb2

bb2:                                              ; preds = %start
  %1 = getelementptr inbounds nuw { ptr, i64 }, ptr %_7, i64 0
  store ptr @alloc_57d70e9d94c65ecfc15225d29a5ed72b, ptr %1, align 8
  %2 = getelementptr inbounds i8, ptr %1, i64 8
  store i64 198, ptr %2, align 8
  store ptr %_7, ptr %_5, align 8
  %3 = getelementptr inbounds i8, ptr %_5, i64 8
  store i64 1, ptr %3, align 8
  %4 = load ptr, ptr @anon.82179aa89ed86be757434e5965125d80.0, align 8
  %5 = load i64, ptr getelementptr inbounds (i8, ptr @anon.82179aa89ed86be757434e5965125d80.0, i64 8), align 8
  %6 = getelementptr inbounds i8, ptr %_5, i64 32
  store ptr %4, ptr %6, align 8
  %7 = getelementptr inbounds i8, ptr %6, i64 8
  store i64 %5, ptr %7, align 8
  %8 = getelementptr inbounds i8, ptr %_5, i64 16
  store ptr inttoptr (i64 8 to ptr), ptr %8, align 8
  %9 = getelementptr inbounds i8, ptr %8, i64 8
  store i64 0, ptr %9, align 8
; call core::panicking::panic_nounwind_fmt
  call void @_ZN4core9panicking18panic_nounwind_fmt17h219ec714db0583beE(ptr align 8 %_5, i1 zeroext false, ptr align 8 %0) #13
  unreachable

bb1:                                              ; preds = %start
  ret void
}

; alloc::raw_vec::RawVecInner<A>::with_capacity_in
; Function Attrs: inlinehint nonlazybind uwtable
define { i64, ptr } @"_ZN5alloc7raw_vec20RawVecInner$LT$A$GT$16with_capacity_in17he4aeaac41456d6dcE"(i64 %capacity, i64 %elem_layout.0, i64 %elem_layout.1, ptr align 8 %0) unnamed_addr #1 {
start:
  %self = alloca [8 x i8], align 8
  %elem_layout = alloca [16 x i8], align 8
  %this = alloca [16 x i8], align 8
  %_4 = alloca [24 x i8], align 8
; call alloc::raw_vec::RawVecInner<A>::try_allocate_in
  call void @"_ZN5alloc7raw_vec20RawVecInner$LT$A$GT$15try_allocate_in17h131528e5e850069eE"(ptr sret([24 x i8]) align 8 %_4, i64 %capacity, i1 zeroext false, i64 %elem_layout.0, i64 %elem_layout.1)
  %_5 = load i64, ptr %_4, align 8
  %1 = trunc nuw i64 %_5 to i1
  br i1 %1, label %bb3, label %bb4

bb3:                                              ; preds = %start
  %2 = getelementptr inbounds i8, ptr %_4, i64 8
  %err.0 = load i64, ptr %2, align 8
  %3 = getelementptr inbounds i8, ptr %2, i64 8
  %err.1 = load i64, ptr %3, align 8
; call alloc::raw_vec::handle_error
  call void @_ZN5alloc7raw_vec12handle_error17h637aa921894539e0E(i64 %err.0, i64 %err.1, ptr align 8 %0) #11
  unreachable

bb4:                                              ; preds = %start
  %4 = getelementptr inbounds i8, ptr %_4, i64 8
  %5 = load i64, ptr %4, align 8
  %6 = getelementptr inbounds i8, ptr %4, i64 8
  %7 = load ptr, ptr %6, align 8
  store i64 %5, ptr %this, align 8
  %8 = getelementptr inbounds i8, ptr %this, i64 8
  store ptr %7, ptr %8, align 8
  store i64 %elem_layout.0, ptr %elem_layout, align 8
  %9 = getelementptr inbounds i8, ptr %elem_layout, i64 8
  store i64 %elem_layout.1, ptr %9, align 8
  %10 = icmp eq i64 %elem_layout.1, 0
  br i1 %10, label %bb6, label %bb7

bb6:                                              ; preds = %bb4
  store i64 -1, ptr %self, align 8
  br label %bb5

bb7:                                              ; preds = %bb4
  %self1 = load i64, ptr %this, align 8
  store i64 %self1, ptr %self, align 8
  br label %bb5

bb5:                                              ; preds = %bb7, %bb6
  %11 = load i64, ptr %self, align 8
  %_13 = sub i64 %11, 0
  %_8 = icmp ugt i64 %capacity, %_13
  %cond = xor i1 %_8, true
  br label %bb8

bb8:                                              ; preds = %bb5
; call core::hint::assert_unchecked::precondition_check
  call void @_ZN4core4hint16assert_unchecked18precondition_check17h979bd9711919bbd4E(i1 zeroext %cond, ptr align 8 %0) #15
  br label %bb9

bb9:                                              ; preds = %bb8
  %_0.0 = load i64, ptr %this, align 8
  %12 = getelementptr inbounds i8, ptr %this, i64 8
  %_0.1 = load ptr, ptr %12, align 8
  %13 = insertvalue { i64, ptr } poison, i64 %_0.0, 0
  %14 = insertvalue { i64, ptr } %13, ptr %_0.1, 1
  ret { i64, ptr } %14

bb2:                                              ; No predecessors!
  unreachable
}

; <T as alloc::slice::<impl [T]>::to_vec_in::ConvertVec>::to_vec
; Function Attrs: inlinehint nonlazybind uwtable
define void @"_ZN87_$LT$T$u20$as$u20$alloc..slice..$LT$impl$u20$$u5b$T$u5d$$GT$..to_vec_in..ConvertVec$GT$6to_vec17h6b2167f36aa25ff3E"(ptr sret([24 x i8]) align 8 %_0, ptr align 1 %s.0, i64 %s.1) unnamed_addr #1 {
start:
  %_21 = alloca [8 x i8], align 8
  %v = alloca [24 x i8], align 8
; call alloc::raw_vec::RawVecInner<A>::with_capacity_in
  %0 = call { i64, ptr } @"_ZN5alloc7raw_vec20RawVecInner$LT$A$GT$16with_capacity_in17he4aeaac41456d6dcE"(i64 %s.1, i64 1, i64 1, ptr align 8 @alloc_ac65f47950d5bfcf522958a117ed0fd0)
  %_11.0 = extractvalue { i64, ptr } %0, 0
  %_11.1 = extractvalue { i64, ptr } %0, 1
  store i64 %_11.0, ptr %v, align 8
  %1 = getelementptr inbounds i8, ptr %v, i64 8
  store ptr %_11.1, ptr %1, align 8
  %2 = getelementptr inbounds i8, ptr %v, i64 16
  store i64 0, ptr %2, align 8
  %3 = getelementptr inbounds i8, ptr %v, i64 8
  %_13 = load ptr, ptr %3, align 8
  br label %bb2

bb2:                                              ; preds = %start
; call core::ptr::copy_nonoverlapping::precondition_check
  call void @_ZN4core3ptr19copy_nonoverlapping18precondition_check17hd1660667b8df50cdE(ptr %s.0, ptr %_13, i64 1, i64 1, i64 %s.1, ptr align 8 @alloc_b72098898095f52adc73e53a589e397d) #15
  %4 = mul i64 %s.1, 1
  call void @llvm.memcpy.p0.p0.i64(ptr align 1 %_13, ptr align 1 %s.0, i64 %4, i1 false)
  br label %bb5

bb5:                                              ; preds = %bb2
  br label %bb10

bb10:                                             ; preds = %bb5
  %self = load i64, ptr %v, align 8
  store i64 %self, ptr %_21, align 8
  br label %bb8

bb9:                                              ; No predecessors!
  store i64 -1, ptr %_21, align 8
  br label %bb8

bb8:                                              ; preds = %bb10, %bb9
  %5 = load i64, ptr %_21, align 8
; call alloc::vec::Vec<T,A>::set_len::precondition_check
  call void @"_ZN5alloc3vec16Vec$LT$T$C$A$GT$7set_len18precondition_check17h36c85061f27ed97bE"(i64 %s.1, i64 %5, ptr align 8 @alloc_2777e455911d18e2781b88be07504468) #15
  br label %bb7

bb7:                                              ; preds = %bb8
  %6 = getelementptr inbounds i8, ptr %v, i64 16
  store i64 %s.1, ptr %6, align 8
  call void @llvm.memcpy.p0.p0.i64(ptr align 8 %_0, ptr align 8 %v, i64 24, i1 false)
  ret void

bb4:                                              ; No predecessors!
  unreachable
}

; probe1::probe
; Function Attrs: nonlazybind uwtable
define void @_ZN6probe15probe17h778e25df45a57be2E() unnamed_addr #3 {
start:
  %_5 = alloca [16 x i8], align 8
  %args = alloca [16 x i8], align 8
  %_3 = alloca [48 x i8], align 8
  %_2 = alloca [24 x i8], align 8
  %_1 = alloca [24 x i8], align 8
; call core::fmt::rt::Argument::new_lower_exp
  call void @_ZN4core3fmt2rt8Argument13new_lower_exp17hdb28f0a75dad489fE(ptr sret([16 x i8]) align 8 %_5, ptr align 8 @alloc_53973d2fe29b4adba8bb7390b5678745)
  %0 = getelementptr inbounds nuw %"core::fmt::rt::Argument<'_>", ptr %args, i64 0
  call void @llvm.memcpy.p0.p0.i64(ptr align 8 %0, ptr align 8 %_5, i64 16, i1 false)
; call core::fmt::rt::<impl core::fmt::Arguments>::new_v1
  call void @"_ZN4core3fmt2rt38_$LT$impl$u20$core..fmt..Arguments$GT$6new_v117hcf32058b4dd71b39E"(ptr sret([48 x i8]) align 8 %_3, ptr align 8 @alloc_eab5d04767146d7d9b93b60d28ef530a, ptr align 8 %args)
; call alloc::fmt::format
  call void @_ZN5alloc3fmt6format17h6ed158a03a5de8a7E(ptr sret([24 x i8]) align 8 %_2, ptr align 8 %_3)
  call void @llvm.memcpy.p0.p0.i64(ptr align 8 %_1, ptr align 8 %_2, i64 24, i1 false)
; call core::ptr::drop_in_place<alloc::string::String>
  call void @"_ZN4core3ptr42drop_in_place$LT$alloc..string..String$GT$17hacbd5aff21b03340E"(ptr align 8 %_1)
  ret void
}

; core::fmt::num::imp::<impl core::fmt::LowerExp for isize>::fmt
; Function Attrs: nonlazybind uwtable
declare zeroext i1 @"_ZN4core3fmt3num3imp55_$LT$impl$u20$core..fmt..LowerExp$u20$for$u20$isize$GT$3fmt17h47b185332d7ed1c2E"(ptr align 8, ptr align 8) unnamed_addr #3

; Function Attrs: nocallback nofree nounwind willreturn memory(argmem: readwrite)
declare void @llvm.memcpy.p0.p0.i64(ptr noalias writeonly captures(none), ptr noalias readonly captures(none), i64, i1 immarg) #4

; Function Attrs: nounwind nonlazybind uwtable
declare i32 @rust_eh_personality(i32, i32, i64, ptr, ptr) unnamed_addr #5

; Function Attrs: nocallback nofree nosync nounwind speculatable willreturn memory(none)
declare i64 @llvm.ctpop.i64(i64) #6

; core::panicking::panic_cannot_unwind
; Function Attrs: cold minsize noinline noreturn nounwind nonlazybind optsize uwtable
declare void @_ZN4core9panicking19panic_cannot_unwind17h48697f70b892f57bE() unnamed_addr #7

; core::panicking::panic_fmt
; Function Attrs: cold noinline noreturn nonlazybind uwtable
declare void @_ZN4core9panicking9panic_fmt17h0023f2c1f721ac1dE(ptr align 8, ptr align 8) unnamed_addr #8

; core::panicking::panic_nounwind_fmt
; Function Attrs: cold noinline noreturn nounwind nonlazybind uwtable
declare void @_ZN4core9panicking18panic_nounwind_fmt17h219ec714db0583beE(ptr align 8, i1 zeroext, ptr align 8) unnamed_addr #9

; <alloc::vec::Vec<T,A> as core::ops::drop::Drop>::drop
; Function Attrs: nonlazybind uwtable
declare void @"_ZN70_$LT$alloc..vec..Vec$LT$T$C$A$GT$$u20$as$u20$core..ops..drop..Drop$GT$4drop17h991b5fea7f8f6006E"(ptr align 8) unnamed_addr #3

; core::panicking::panic_in_cleanup
; Function Attrs: cold minsize noinline noreturn nounwind nonlazybind optsize uwtable
declare void @_ZN4core9panicking16panic_in_cleanup17h7a69428c8f5fb9a4E() unnamed_addr #7

; <alloc::raw_vec::RawVec<T,A> as core::ops::drop::Drop>::drop
; Function Attrs: nonlazybind uwtable
declare void @"_ZN77_$LT$alloc..raw_vec..RawVec$LT$T$C$A$GT$$u20$as$u20$core..ops..drop..Drop$GT$4drop17h902845a8e4723e86E"(ptr align 8) unnamed_addr #3

; Function Attrs: nocallback nofree nosync nounwind speculatable willreturn memory(none)
declare { i64, i1 } @llvm.umul.with.overflow.i64(i64, i64) #6

; core::panicking::panic_nounwind
; Function Attrs: cold noinline noreturn nounwind nonlazybind uwtable
declare void @_ZN4core9panicking14panic_nounwind17h47fee6258c2f7571E(ptr align 1, i64) unnamed_addr #9

; alloc::fmt::format::format_inner
; Function Attrs: nonlazybind uwtable
declare void @_ZN5alloc3fmt6format12format_inner17h52d1570151877af0E(ptr sret([24 x i8]) align 8, ptr align 8) unnamed_addr #3

; alloc::raw_vec::RawVecInner<A>::try_allocate_in
; Function Attrs: nonlazybind uwtable
declare void @"_ZN5alloc7raw_vec20RawVecInner$LT$A$GT$15try_allocate_in17h131528e5e850069eE"(ptr sret([24 x i8]) align 8, i64, i1 zeroext, i64, i64) unnamed_addr #3

; alloc::raw_vec::handle_error
; Function Attrs: cold minsize noreturn nonlazybind optsize uwtable
declare void @_ZN5alloc7raw_vec12handle_error17h637aa921894539e0E(i64, i64, ptr align 8) unnamed_addr #10

attributes #0 = { cold nounwind nonlazybind uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #1 = { inlinehint nonlazybind uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #2 = { inlinehint nounwind nonlazybind uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #3 = { nonlazybind uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #4 = { nocallback nofree nounwind willreturn memory(argmem: readwrite) }
attributes #5 = { nounwind nonlazybind uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #6 = { nocallback nofree nosync nounwind speculatable willreturn memory(none) }
attributes #7 = { cold minsize noinline noreturn nounwind nonlazybind optsize uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #8 = { cold noinline noreturn nonlazybind uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #9 = { cold noinline noreturn nounwind nonlazybind uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #10 = { cold minsize noreturn nonlazybind optsize uwtable "probe-stack"="inline-asm" "target-cpu"="x86-64" }
attributes #11 = { noreturn }
attributes #12 = { cold noreturn nounwind }
attributes #13 = { noreturn nounwind }
attributes #14 = { cold }
attributes #15 = { nounwind }

!llvm.module.flags = !{!0, !1}
!llvm.ident = !{!2}

!0 = !{i32 8, !"PIC Level", i32 2}
!1 = !{i32 2, !"RtLibUseGOT", i32 1}
!2 = !{!"rustc version 1.91.0-nightly (8e62bfd31 2025-08-12)"}
