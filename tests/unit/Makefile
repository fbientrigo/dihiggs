# tests/unit/Makefile

# Compilador y Flags (Copiamos lo básico de tu Makefile principal)
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O0 -g
# Importante: Incluir la carpeta src del proyecto principal
INCLUDES = -I../../dihiggs/src

# Librerías necesarias (Matemáticas es vital para ParamUtils)
LDFLAGS = -lm 

# Si en el futuro ParamUtils usa GSL, descomenta esto:
# LDFLAGS += -lgsl -lgslcblas

# Archivos
TEST_SRCS = test_main.cpp test_ParamUtils.cpp
# El código fuente que estamos probando (fuera de la carpeta de tests)
IMPL_SRCS = ../../dihiggs/src/ParamUtils.cpp

# Binario de salida
TARGET = run_tests

all: $(TARGET)

$(TARGET): $(TEST_SRCS) $(IMPL_SRCS)
	@echo "🧪 Compilando Tests Unitarios..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@ $(LDFLAGS)
	@echo "✅ Compilación exitosa. Ejecuta ./$(TARGET)"

clean:
	rm -f $(TARGET)