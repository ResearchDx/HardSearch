################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AlignmentNode.cpp \
../src/BedUtils.cpp \
../src/BreakPointNode.cpp \
../src/HashEntry.cpp \
../src/HashMap.cpp \
../src/MurmurHash3.cpp \
../src/RapidXMLWrapper.cpp \
../src/ReadSequence.cpp \
../src/SamParser.cpp \
../src/SoftClipRead.cpp \
../src/Utils.cpp \
../src/main.cpp \
../src/blast_helper.cpp

OBJS += \
./src/AlignmentNode.o \
./src/BedUtils.o \
./src/BreakPointNode.o \
./src/HashEntry.o \
./src/HashMap.o \
./src/MurmurHash3.o \
./src/RapidXMLWrapper.o \
./src/ReadSequence.o \
./src/SamParser.o \
./src/SoftClipRead.o \
./src/Utils.o \
./src/main.o \

BH_OBJS += \
./src/blast_helper.o


CPP_DEPS += \
./src/AlignmentNode.d \
./src/BedUtils.d \
./src/BreakPointNode.d \
./src/HashEntry.d \
./src/HashMap.d \
./src/MurmurHash3.d \
./src/RapidXMLWrapper.d \
./src/ReadSequence.d \
./src/SamParser.d \
./src/SoftClipRead.d \
./src/Utils.d \
./src/main.d


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../lib/samtools-0.1.18 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


