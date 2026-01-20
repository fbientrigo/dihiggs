#include "catch.hpp"
#include "../../dihiggs/src/ParamUtils.hpp"
#include <cmath>

// Definimos un margen de error pequeño para comparar doubles
#define DOUBLE_TOLERANCE 1e-9

TEST_CASE("Utilidades de Bucle (stepsCount)", "[core][utils]") {
    // Caso: Rango exacto (0 a 10 con pasos de 1) -> 11 puntos (0,1,..,10)
    REQUIRE(stepsCount(0.0, 10.0, 1.0) == 11);

    // Caso: Rango decimal
    REQUIRE(stepsCount(0.0, 1.0, 0.1) == 11);

    // Caso: Un solo punto (min == max)
    REQUIRE(stepsCount(5.0, 5.0, 1.0) == 1);
    
    // Caso: Paso más grande que el rango (debe dar al menos el punto inicial y final si ajustamos la logica, 
    // pero tu implementacion usa floor, asi que verifiquemos su comportamiento actual)
    // Tu lógica: floor((10-0)/20) + 1 = 0 + 1 = 1. Correcto.
    REQUIRE(stepsCount(0.0, 10.0, 20.0) == 1);
}

TEST_CASE("Física: Chequeo de Positividad (Stability)", "[physics][core]") {
    
    SECTION("Caso Válido (Todos positivos)") {
        // Lambda1..5 positivos cumplen trivialmente
        REQUIRE(check_positivity(1.0, 1.0, 1.0, 1.0, 0.0) == true);
    }

    SECTION("Caso Fallido: Lambda1 negativo") {
        REQUIRE(check_positivity(-1.0, 1.0, 1.0, 1.0, 0.0) == false);
    }

    SECTION("Caso Fallido: Lambda2 negativo") {
        REQUIRE(check_positivity(1.0, -1.0, 1.0, 1.0, 0.0) == false);
    }

    SECTION("Caso Límite: Lambda3 viola condición") {
        // Condición: lambda3 > -sqrt(l1*l2)
        // sqrt(1*1) = 1. Si lambda3 es -2, falla.
        REQUIRE(check_positivity(1.0, 1.0, -2.0, 0.0, 0.0) == false);
    }
}

TEST_CASE("Física: Perturbatividad (check_lambda)", "[physics]") {
    // Límite es 4*PI^2 ~ 157.9 o |lambda| < 4*PI (~12.56)
    // Tu funcion chequea: l*l < 16*pi*pi
    
    //double limit = 4.0 * 3.14159265359; 

    SECTION("Dentro del límite") {
        REQUIRE(check_lambda(1.0) == true);
        REQUIRE(check_lambda(10.0) == true);
    }

    SECTION("Fuera del límite") {
        REQUIRE(check_lambda(20.0) == false); // 20 > 12.56
        REQUIRE(check_lambda(-20.0) == false);
    }
}

TEST_CASE("Cálculos: Lambda 1", "[physics][math]") {
    // Probamos calc_lambda1 con valores conocidos o neutros para ver que no explota
    // mh=125, mH=300, m12_2=100^2, sa=0, ca=1, sb=1, cb=0... cuidado con dividir por cero
    // Tu funcion divide por (VEV*VEV * cb*cb). Si cb es 0, explota.
    
    double mh = 125.0;
    double mH = 300.0;
    double m12_2 = 10000.0;
    double sa = 0.1;
    double ca = std::sqrt(1 - sa*sa);
    double tan_beta = 10.0;
    
    // Recalculamos cb y sb basados en tan_beta para consistencia
    double cb = 1.0 / std::sqrt(1 + tan_beta*tan_beta); 
    double sb = tan_beta * cb;

    double lam6 = 0.0;
    double lam7 = 0.0;
    double VEV = 246.0;

    double result = calc_lambda1(mh, mH, m12_2, sa, ca, sb, cb, tan_beta, lam6, lam7, VEV);

    // Solo verificamos que sea un número finito (no NaN ni Inf)
    REQUIRE_FALSE(std::isnan(result));
    REQUIRE_FALSE(std::isinf(result));
}