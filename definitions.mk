# written by Chihun Sung & Shyamini Kularathna

WORK_DIR = /home/chihun/MPM/u-p-rigid2

# compiler/linker
C_PP = g++

# compilation/linking flags
C_FLAGS = -Wall -pedantic -fopenmp -pthread 
L_FLAGS = -lm -lpthread -lgomp 

# library paths and names
LIB_DIR = $(WORK_DIR)/libs
BOOST = $(LIB_DIR)/boost_1_73_0
EIGEN = /usr/local/include/eigen3

# include directory for mpm header files
INCLUDE_DIR = $(WORK_DIR)/includes 
MATERIAL = $(WORK_DIR)/includes/material
SRC = $(WORK_DIR)/includes/src
MISC = $(SRC)/misc
MESH = $(SRC)/mesh
PARTICLE = $(SRC)/particle
MATRIX = $(SRC)/matrix
SOLVER = $(SRC)/solvers
