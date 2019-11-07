################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Analysis.cpp \
../src/CrankNicolson.cpp \
../src/FirstOrderExplicit.cpp \
../src/Main.cpp \
../src/Mesh.cpp \
../src/Point3D.cpp 

OBJS += \
./src/Analysis.o \
./src/CrankNicolson.o \
./src/FirstOrderExplicit.o \
./src/Main.o \
./src/Mesh.o \
./src/Point3D.o 

CPP_DEPS += \
./src/Analysis.d \
./src/CrankNicolson.d \
./src/FirstOrderExplicit.d \
./src/Main.d \
./src/Mesh.d \
./src/Point3D.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I"/home/aishu/Documents/Scalar-Transport-Equation-Solver/include" -O0 -g3 -Wall -Werror -c -fmessage-length=0 -fsplit-stack -fPIC -pthread -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


