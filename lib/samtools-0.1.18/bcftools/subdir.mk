################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../lib/samtools-0.1.18/bcftools/bcf.o \
../lib/samtools-0.1.18/bcftools/bcf2qcall.o \
../lib/samtools-0.1.18/bcftools/bcfutils.o \
../lib/samtools-0.1.18/bcftools/call1.o \
../lib/samtools-0.1.18/bcftools/em.o \
../lib/samtools-0.1.18/bcftools/fet.o \
../lib/samtools-0.1.18/bcftools/index.o \
../lib/samtools-0.1.18/bcftools/kfunc.o \
../lib/samtools-0.1.18/bcftools/kmin.o \
../lib/samtools-0.1.18/bcftools/main.o \
../lib/samtools-0.1.18/bcftools/mut.o \
../lib/samtools-0.1.18/bcftools/prob1.o \
../lib/samtools-0.1.18/bcftools/vcf.o 

C_SRCS += \
../lib/samtools-0.1.18/bcftools/bcf.c \
../lib/samtools-0.1.18/bcftools/bcf2qcall.c \
../lib/samtools-0.1.18/bcftools/bcfutils.c \
../lib/samtools-0.1.18/bcftools/call1.c \
../lib/samtools-0.1.18/bcftools/em.c \
../lib/samtools-0.1.18/bcftools/fet.c \
../lib/samtools-0.1.18/bcftools/index.c \
../lib/samtools-0.1.18/bcftools/kfunc.c \
../lib/samtools-0.1.18/bcftools/kmin.c \
../lib/samtools-0.1.18/bcftools/main.c \
../lib/samtools-0.1.18/bcftools/mut.c \
../lib/samtools-0.1.18/bcftools/prob1.c \
../lib/samtools-0.1.18/bcftools/vcf.c 

OBJS += \
./lib/samtools-0.1.18/bcftools/bcf.o \
./lib/samtools-0.1.18/bcftools/bcf2qcall.o \
./lib/samtools-0.1.18/bcftools/bcfutils.o \
./lib/samtools-0.1.18/bcftools/call1.o \
./lib/samtools-0.1.18/bcftools/em.o \
./lib/samtools-0.1.18/bcftools/fet.o \
./lib/samtools-0.1.18/bcftools/index.o \
./lib/samtools-0.1.18/bcftools/kfunc.o \
./lib/samtools-0.1.18/bcftools/kmin.o \
./lib/samtools-0.1.18/bcftools/main.o \
./lib/samtools-0.1.18/bcftools/mut.o \
./lib/samtools-0.1.18/bcftools/prob1.o \
./lib/samtools-0.1.18/bcftools/vcf.o 

C_DEPS += \
./lib/samtools-0.1.18/bcftools/bcf.d \
./lib/samtools-0.1.18/bcftools/bcf2qcall.d \
./lib/samtools-0.1.18/bcftools/bcfutils.d \
./lib/samtools-0.1.18/bcftools/call1.d \
./lib/samtools-0.1.18/bcftools/em.d \
./lib/samtools-0.1.18/bcftools/fet.d \
./lib/samtools-0.1.18/bcftools/index.d \
./lib/samtools-0.1.18/bcftools/kfunc.d \
./lib/samtools-0.1.18/bcftools/kmin.d \
./lib/samtools-0.1.18/bcftools/main.d \
./lib/samtools-0.1.18/bcftools/mut.d \
./lib/samtools-0.1.18/bcftools/prob1.d \
./lib/samtools-0.1.18/bcftools/vcf.d 


# Each subdirectory must supply rules for building sources it contributes
lib/samtools-0.1.18/bcftools/%.o: ../lib/samtools-0.1.18/bcftools/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


