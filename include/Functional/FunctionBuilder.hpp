#pragma once 

#ifndef DMATH_EUQATION_HPP
#define DMATH_EUQATION_HPP

#include"mathFunction.hpp"
#include"../Utils/dataTypes.hpp"
#include"../Complex/complexFunctions.hpp"
#include"vectorfunction.hpp"


NAMESPACESTART

/* this file contains a builder pattern for specific equations
 * like: quadratic equations, linear equations, etc.
 * This header is part of the source file: @filr/mathFunction.hpp
*/

#pragma region MathFunctionBuilders


Dmath::SingleVarFunction sinBuilder(Dmath::Scalar a, Dmath::Scalar b, Dmath::Scalar c) {
    std::string functionString =
        std::to_string(a) + " * sin(" +
        std::to_string(b) + " * x";

    if(c > 0) {
        functionString += " + " + std::to_string(c);
    }
    else if(c < 0) {
        functionString += " - " + std::to_string(std::abs(c));
    }

    functionString += ")";

    return Dmath::SingleVarFunction(
        [a, b, c](double x) -> double {
            return a * std::sin(b * x + c);
        },
        functionString
    );
}

Dmath::SingleVarFunction cosBuilder(Dmath::Scalar a, Dmath::Scalar b, Dmath::Scalar c) {
    std::string functionString =
        std::to_string(a) + " * cos(" +
        std::to_string(b) + " * x";

    if(c > 0) {
        functionString += " + " + std::to_string(c);
    }
    else if(c < 0) {
        functionString += " - " + std::to_string(std::abs(c));
    }

    functionString += ")";

    return Dmath::SingleVarFunction(
        [a, b, c](double x) -> double {
            return a * std::cos(b * x + c);
        },
        functionString
    );
}

Dmath::SingleVarFunction tanBuilder(Dmath::Scalar a, Dmath::Scalar b, Dmath::Scalar c) {
    std::string functionString =
        std::to_string(a) + " * tan(" +
        std::to_string(b) + " * x";

    if(c > 0) {
        functionString += " + " + std::to_string(c);
    }
    else if(c < 0) {
        functionString += " - " + std::to_string(std::abs(c));
    }

    functionString += ")";

    return Dmath::SingleVarFunction(
        [a, b, c](double x) -> double {
            return a * std::tan(b * x + c);
        },
        functionString
    );
}

Dmath::SingleVarFunction cotBuilder(Dmath::Scalar a, Dmath::Scalar b, Dmath::Scalar c) {
    std::string functionString =
        std::to_string(a) + " * cot(" +
        std::to_string(b) + " * x";

    if(c > 0) {
        functionString += " + " + std::to_string(c);
    }
    else if(c < 0) {
        functionString += " - " + std::to_string(std::abs(c));
    }

    functionString += ")";

    return Dmath::SingleVarFunction(
        [a, b, c](double x) -> double {
            return a * Dmath::cot(b * x + c);
        },
        functionString
    );
}


//Linear equation builder
// This function returns a lambda that represents a linear equation of the form y = ax + b
Dmath::SingleVarFunction linerarEquation(Dmath::Scalar a, Dmath::Scalar b) {
    std::string functionString = std::to_string(a) + " * x +" + std::to_string(b);
    return Dmath::SingleVarFunction([a, b](double x) -> double {
        return a * x + b;
    },functionString);
}

//Quadratic equation builder in its vertex form: f(x) = a*(x+b)² +c
Dmath::SingleVarFunction quadraticEquationVertexForm(Dmath::Scalar a, Dmath::Scalar b, Dmath::Scalar c) {
    std::string functionString = std::to_string(a) + " * (x +" + std::to_string(b) + ")^2 + " + std::to_string(c);

    return Dmath::SingleVarFunction(  [a, b, c](double x) -> double {
        return a * std::pow((x+b), 2)  + c;
    }, functionString);
}

Dmath::SingleVarFunction quadraticEquationNormalForm(Dmath::Scalar a, Dmath::Scalar b, Dmath::Scalar c) {
    std::string functionString = std::to_string(a) + " * x^2 +" + std::to_string(b) + "x + " + std::to_string(c);
    
    return Dmath::SingleVarFunction( [a, b, c](double x) -> double {
        return a * std::pow(x, 2) + b * x  + c;
    },functionString);
}

//cubic equation builder
// y = ax³ + bx² + cx + d
Dmath::SingleVarFunction cubicEquation(Dmath::Scalar a, Dmath::Scalar b, Dmath::Scalar c, Dmath::Scalar d) {

    std::string functionString = std::to_string(a) + " * x^3 +" + std::to_string(b) + "* x^2 +" + std::to_string(b) + "* x^2" ;
    return Dmath::SingleVarFunction( [a, b, c, d](double x) -> double {
        return a * x * x * x + b * x * x + c * x + d;
    
    },functionString);
    
}

//Exponential equation builder
// y = a * b^x
Dmath::SingleVarFunction exponentialEquation(Dmath::Scalar a, Dmath::Scalar b) {
    std::string functionString = std::to_string(a) + " * "+  std::to_string(b) +"^x";
    return  Dmath::SingleVarFunction ([a, b](double x) -> double {
        return a * std::pow(b, x);
    }, functionString);
}

//Logarithmic equation builder
// y = a * log_b(x)
Dmath::SingleVarFunction logarithmicEquation(Dmath::Scalar a, Dmath::Scalar b) {


    std::string functionString = std::to_string(a) + " * log_b(x)";

    return [a, b](double x) -> double {
        Dmath::Scalar logValue = a * Dmath::logarithm(b,x);
        if (logValue == 0) {
            throw std::invalid_argument("Logarithm of zero is undefined.");
            return 0;
        }
        return logValue;
    };
}



//jacobiDeterminantSpherecoordinates
// This function returns a lambda that calculates the Jacobian determinant for spherical coordinates
Dmath::DoubleVarFunction jacobiDeterminantSphere(){
    return [](double r, double theta) ->double {
        // Jacobian determinant for spherical coordinates: r^2 * sin(theta)
        return r * r * std::sin(theta);
    };
}

// Jacobi determinant for cylindrical coordinates
Dmath::DoubleVarFunction jacobiDeterminantCylindrical() {
    return [](double r, double theta) -> double {
        // Jacobian determinant for cylindrical coordinates: r
        return r;
    };
}

#pragma endregion //MathFunctionBuilders


#pragma region VectorFunctionBuilders

Dmath::SingleVectorFunction circleEquation(Dmath::Scalar radius) {
    return Dmath::SingleVectorFunction(
        Dmath::SingleVarFunction(
            [radius](double t) -> double {
                return radius * std::cos(t); 
            },std::to_string(radius) + "cos(t)"), // x(t)

        Dmath::SingleVarFunction(
            [radius](double t) -> double {
                return radius * std::sin(t); 
            },std::to_string(radius) + "sin(t)")
    );
}


Dmath::SingleVectorFunction unitCircle() {
    return Dmath::SingleVectorFunction(
        [](double t) -> double { return std::cos(t); }, // x(t)
        [](double t) -> double { return std::sin(t); }  // y(t)
    );
}

Dmath::SingleVectorFunction elipseEquation(Dmath::Scalar a, Dmath::Scalar b) {
    return Dmath::SingleVectorFunction(
        Dmath::SingleVarFunction([a](double t) -> double { return a * std::cos(t); }, std::to_string(a) + "* cos(t)"),
        
        Dmath::SingleVarFunction([b](double t) -> double { return b * std::sin(t); }, std::to_string(b) + "* sin(t)")
    );
}


Dmath::SingleVectorFunction helixEquation(Dmath::Scalar radius, Dmath::Scalar height) {
    return Dmath::SingleVectorFunction(
        [radius](double t) -> double { return radius * std::cos(t); }, // x(t)
        [radius](double t) -> double { return radius * std::sin(t); }, // y(t)
        [height](double t) -> double { return height * t; }            // z(t)
    );
}
// x(t) = a * cosh(t/a)
// y(t) = t 
Dmath::SingleVectorFunction ropeCurve(Dmath::Scalar a){
    return Dmath::SingleVectorFunction(
        [a](double t) -> double { return a * std::cosh(t / a); }, // x(t)
        [](double t) -> double { return t; }                      // y(t)
    );
}


//Wurfparabel
Dmath::SingleVectorFunction parabolicTrajectory(Dmath::Scalar v, Dmath::Scalar angle, double height){
    Dmath::SingleVarFunction xOfT = [v,angle,height](double t){
        return v * std::cos(angle) * t; // x(t)
    };
    Dmath::SingleVarFunction yOfT = [v,angle,height](double t){
        return   height + v * std::sin(angle) * t - 0.5* 9.81 * t   * t ; // y(t)
    };
    return Dmath::SingleVectorFunction(xOfT, yOfT);
}

Dmath::DoubleVectorFunction sphereEquation(Dmath::Scalar radius) {
    return Dmath::DoubleVectorFunction(
        [radius](double u, double v) -> double { return radius * std::sin(u) * std::cos(v); }, // x(u, v)
        [radius](double u, double v) -> double { return radius * std::sin(u) * std::sin(v); }, // y(u, v)
        [radius](double u, double v) -> double { return radius * std::cos(u); }                // z(u, v)
    );
}


#pragma endregion //VectorFunctionBuilders






Dmath::SingleVarFunction parseTrigFunction(std::string funcName, const std::string& input) {
    std::string str = input;
    std::string efn = funcName + "(";

    str = Dmath::StringHelper().removeWhitespace(str);

    if(str.size() < 6) {
        throw std::invalid_argument(
            "Build Error: Functional/FunctionBuilder.hpp parseSinFunction(const std::string& input):\n Invalid sine function: " + input
        );
    }

    if(str.rfind(efn, 0) == 0) {
        // sin(...)
        str = "1*" + str;
    }

    size_t sinPos = str.find(efn);

    if(sinPos == std::string::npos) {
        throw std::invalid_argument(
            "Expected f(...): " + input
        );
    }

    std::string coefficientString = str.substr(0, sinPos);

    if(coefficientString.size() < 2 ||
       coefficientString.back() != '*') {
        throw std::invalid_argument(
            "Invalid sinus coefficient: " + input
        );
    }

    coefficientString.pop_back();

    Dmath::Scalar a = 1.0;

    if(coefficientString == "-") {
        a = -1.0;
    }
    else if(coefficientString != "+") {
        try {
            a = std::stod(coefficientString);
        }
        catch(...) {
            throw std::invalid_argument(
                "Invalid sinus coefficient: " + coefficientString
            );
        }
    }

    size_t bracketStart = str.find('(', sinPos);
    size_t bracketEnd = str.rfind(')');

    if(bracketStart == std::string::npos ||
       bracketEnd == std::string::npos ||
       bracketEnd <= bracketStart) {
        throw std::invalid_argument(
            "Invalid sin brackets: " + input
        );
    }

    std::string inside =
        str.substr(
            bracketStart + 1,
            bracketEnd - bracketStart - 1
        );

    Dmath::Scalar b = 1.0;
    Dmath::Scalar c = 0.0;

    if(inside == "x") {
        if(funcName == "sin"){
            return Dmath::sinBuilder(a, b, c);
        }
        if(funcName == "cos"){
            return Dmath::cosBuilder(a, b, c);
        }
    }

    size_t xPos = inside.find('x');

    if(xPos == std::string::npos) {
        throw std::invalid_argument(
            "Expected x inside sin(): " + input
        );
    }

    std::string beforeX = inside.substr(0, xPos);
    std::string afterX = inside.substr(xPos + 1);

    // 2*x
    if(beforeX.size() >= 2 && beforeX.back() == '*') {
        beforeX.pop_back();

        if(!beforeX.empty()) {
            b = std::stod(beforeX);
        }
    }

    // x + c / x - c
    if(!afterX.empty()) {
        if(afterX[0] == '+') {
            c = std::stod(afterX.substr(1));
        }
        else if(afterX[0] == '-') {
            c = -std::stod(afterX.substr(1));
        }
        else {
            throw std::invalid_argument(
                "Build Error: .../Functional/FunctionBuilder.hpp Dmath::SingleVarFunction parseTrigFunction(std::string funcName,const std::string& input) : \n\t "
                "Invalid function offset: " + afterX
            );
        }
    }

    if(funcName == "sin") {
        return Dmath::sinBuilder(a, b, c);
    }
    else if(funcName == "cos") {
        return Dmath::cosBuilder(a, b, c);
    }
    else if(funcName == "tan") {
        return Dmath::tanBuilder(a, b, c);
    }
    else if(funcName == "cot") {
        return Dmath::cotBuilder(a, b, c);
    }

    throw std::invalid_argument(
        "Unknown trigonometric function: " + funcName
    );
}

NAMESPACEEND

#endif //DMATH_EUQATION_HPP