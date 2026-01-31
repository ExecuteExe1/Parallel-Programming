	.file	"ArgumentsToThreads.c"
	.text
	.globl	primes
	.data
	.align 32
	.type	primes, @object
	.size	primes, 40
primes:
	.long	2
	.long	3
	.long	5
	.long	7
	.long	11
	.long	13
	.long	17
	.long	19
	.long	23
	.long	29
	.section	.rodata
.LC0:
	.string	"%d"
	.text
	.globl	routine
	.type	routine, @function
routine:
.LFB6:
	.cfi_startproc
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$32, %rsp
	movq	%rdi, -24(%rbp)
	movq	-24(%rbp), %rax
	movl	(%rax), %eax
	movl	%eax, -4(%rbp)
	movl	-4(%rbp), %eax
	cltq
	leaq	0(,%rax,4), %rdx
	leaq	primes(%rip), %rax
	movl	(%rdx,%rax), %eax
	movl	%eax, %esi
	leaq	.LC0(%rip), %rax
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf@PLT
	nop
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE6:
	.size	routine, .-routine
	.section	.rodata
.LC1:
	.string	"Failed to create thread"
.LC2:
	.string	"Failed to join thread"
	.text
	.globl	main
	.type	main, @function
main:
.LFB7:
	.cfi_startproc
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	addq	$-128, %rsp
	movl	%edi, -116(%rbp)
	movq	%rsi, -128(%rbp)
	movq	%fs:40, %rax
	movq	%rax, -8(%rbp)
	xorl	%eax, %eax
	movl	$0, -104(%rbp)
	jmp	.L3
.L5:
	movl	-104(%rbp), %edx
	leaq	-96(%rbp), %rax
	movslq	%edx, %rdx
	salq	$3, %rdx
	leaq	(%rax,%rdx), %rdi
	leaq	-104(%rbp), %rax
	movq	%rax, %rcx
	leaq	routine(%rip), %rax
	movq	%rax, %rdx
	movl	$0, %esi
	call	pthread_create@PLT
	testl	%eax, %eax
	je	.L4
	leaq	.LC1(%rip), %rax
	movq	%rax, %rdi
	call	perror@PLT
.L4:
	movl	-104(%rbp), %eax
	addl	$1, %eax
	movl	%eax, -104(%rbp)
.L3:
	movl	-104(%rbp), %eax
	cmpl	$9, %eax
	jle	.L5
	movl	$0, -100(%rbp)
	jmp	.L6
.L8:
	movl	-100(%rbp), %eax
	cltq
	movq	-96(%rbp,%rax,8), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	pthread_join@PLT
	testl	%eax, %eax
	je	.L7
	leaq	.LC2(%rip), %rax
	movq	%rax, %rdi
	call	perror@PLT
.L7:
	addl	$1, -100(%rbp)
.L6:
	cmpl	$9, -100(%rbp)
	jle	.L8
	movl	$0, %eax
	movq	-8(%rbp), %rdx
	subq	%fs:40, %rdx
	je	.L10
	call	__stack_chk_fail@PLT
.L10:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE7:
	.size	main, .-main
	.ident	"GCC: (Ubuntu 13.2.0-23ubuntu4) 13.2.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
