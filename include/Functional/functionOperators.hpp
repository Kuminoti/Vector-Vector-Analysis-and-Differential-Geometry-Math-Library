#pragma once 

#ifndef DMATH_FUNCTION_OPERATORS_HPP
#define DMATH_FUNCTION_OPERATORS_HPP

#include"../Text/chars.hpp"
#include"vectorfunction.hpp"
#include"../Types/StringHelper.hpp"


NAMESPACESTART

#pragma region SingleVarOperators

Dmath::SingleVarFunction operator*(const char* op,const Dmath::SingleVarFunction& func) {
    Dmath::StringHelper strH;
    std::string operation(op);

    // Numeric operator: "3" * f
    if(strH.isNumeric(operation)){
        return func * std::stod(operation);
    }

    // Derivative
    if(operation == Dmath::SpecialChars.getFromKey("nabla") ||
       operation == Dmath::SpecialChars.getFromKey("derivative") ||
       operation == Dmath::SpecialChars.getFromKey("partial")){

        Dmath::SingleVarFunction f = func.getDerivative();

        f.setFuncData(
            "d/dx(" + func.getFunctionData() + ")"
        );

        return f;
    }

    // Integral
    if(operation == Dmath::SpecialChars.getFromKey("integral")){

        Dmath::SingleVarFunction f = func.getAntiDerivative();

        f.setFuncData(
            Dmath::SpecialChars.getFromKey("integral") +
            "(" + func.getFunctionData() + ")"
        );

        return f;
    }

    // sqrt
    if(operation == Dmath::SpecialChars.getFromKey("sqrt")){

        return Dmath::SingleVarFunction(
            [func](Dmath::Scalar x) -> Dmath::Scalar {
                return std::sqrt(func(x));
            },
            Dmath::SpecialChars.getFromKey("sqrt") +
            "(" + func.getFunctionData() + ")"
        );
    }

    // cbrt
    if(operation == Dmath::SpecialChars.getFromKey("cbrt")){

        return Dmath::SingleVarFunction(
            [func](Dmath::Scalar x) -> Dmath::Scalar {
                return std::cbrt(func(x));
            },
            Dmath::SpecialChars.getFromKey("cbrt") +
            "(" + func.getFunctionData() + ")"
        );
    }

    // fourth root
    if(operation == Dmath::SpecialChars.getFromKey("fourthRoot")){

        return Dmath::SingleVarFunction(
            [func](Dmath::Scalar x) -> Dmath::Scalar {
                return std::pow(func(x), 1.0 / 4.0);
            },
            Dmath::SpecialChars.getFromKey("fourthRoot") +
            "(" + func.getFunctionData() + ")"
        );
    }

    throw std::invalid_argument(
        "Unknown function operator: " + operation
    );
}



template<typename Dmath_type>
Dmath::SingleVectorFunction operator*(const Dmath_type& data, const Dmath::SingleVectorFunction& func) {

    
    if constexpr(std::is_same<Dmath_type, std::string>::value) {

        if(data == Dmath::SpecialChars.getFromKey("nabla")) {

            return Dmath::SingleVectorFunction(
                func.getXFunc().getDerivative(),
                func.getYFunc().getDerivative(),
                func.getZFunc().getDerivative()
            );
        }
    }


    throw std::invalid_argument(
        "Invalid operator for SingleVectorFunction"
    );
}




    #pragma endregion


#pragma region DoubleVar OP














#pragma endregion //DoubleVar OP

NAMESPACEEND

#endif //DMATH_FUNCTION_OPERATORS_HPP