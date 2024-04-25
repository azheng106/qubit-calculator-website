import math
import re
import cmath
from enum import Enum

import numpy as np


class Precedence(Enum):
    # Higher values indicate higher precedence
    FUNCTION = 4
    UNARY_MINUS = 3
    POWER = 2
    MULTIPLY = 1
    ADD = 0
    OPEN_PARENTHESIS = -99  # When comparing any operator to '(', '(' is always lower precedence


precedences = {
    'sin': Precedence.FUNCTION.value,
    'cos': Precedence.FUNCTION.value,
    'tan': Precedence.FUNCTION.value,
    'sqrt': Precedence.FUNCTION.value,
    'unary_minus': Precedence.UNARY_MINUS.value,
    '^': Precedence.POWER.value,
    '*': Precedence.MULTIPLY.value,
    '/': Precedence.MULTIPLY.value,
    '+': Precedence.ADD.value,
    '-': Precedence.ADD.value,
    '(': Precedence.OPEN_PARENTHESIS.value,
}


def stack_has_higher_precedence(stack, current_token) -> bool:
    """
    Return True if the operator on top of the stack has >= precedence than the current token, False otherwise
    """
    stack_precedence = precedences.get(stack[-1], float('inf'))
    current_precedence = precedences.get(current_token, float('inf'))
    return stack_precedence >= current_precedence


def is_negative_sign(tokens, index) -> bool:
    """
    Determine if the "-" is a negative number (True) or minus sign (False)
    :param index: index of the "-"
    :param tokens: list of tokens
    """
    if index == 0 or tokens[index - 1] in '*/+-^(':
        return True
    return False


def is_float(token) -> bool:
    """
    Check if the input token is a float/int because isnumeric() only works for ints
    """
    if token.replace('.', '').isnumeric():
        return True
    return False


def infix_to_postfix(infix) -> list:
    """
    Use Shunting Yard Algorithm to translate mathematical infix to postfix
    https://www.youtube.com/watch?v=Wz85Hiwi5MY
    :param infix: Input mathematical expression
    :return: Postfix
    """
    infix = infix.replace(' ', '')  # Remove spaces
    tokens = re.findall(r'[\d.]+|[+\-*/()^]|sqrt|i|j|sin|cos|tan|pi|e',
                        infix)  # infix = '1+1' => tokens = ['1', '+', '1']

    reconstructed = ''.join(tokens)
    if reconstructed != infix:  # Without this, an infix of "iabcd" will be tokenized into ['i'] and seen as valid
        raise ValueError("Invalid characters detected in infix expression")

    if tokens.count('(') - 1 == tokens.count(')'):  # If user forgot to add a closing parenthesis
        tokens.append(')')

    for i in range(len(tokens) - 1):  # If user inputs '3i', make sure postfix evaluator knows to multiply 3 and i
        curr_token = tokens[i]
        next_token = tokens[i + 1]
        if is_float(curr_token) and next_token in ['i', 'j', 'pi', 'sqrt', 'sin', 'cos',
                                                   'tan', 'e']:  # ['2', 'i'] => ['2', '*', 'i']
            tokens.insert(i + 1, '*')

    stack = []  # Operator stack
    queue = []  # Output queue

    for index, token in enumerate(tokens):
        if is_float(token) or token in ['i', 'j', 'pi', 'e']:  # Numbers are automatically pushed into queue
            queue.append(token)
        elif token == '-':
            if is_negative_sign(tokens, index):  # Treat negative sign differently
                stack.append('unary_minus')
            else:  # Treat like a normal operator
                if stack and stack_has_higher_precedence(stack, token):
                    queue.append(stack.pop())
                stack.append(token)
        elif token == ')':  # Pop stack into queue until '(' is found, and then discard both parenthesis
            while stack and stack[-1] != '(':
                queue.append(stack.pop())
            stack.pop()  # Pop the '('
        elif token == '(':  # Same as "else" case
            stack.append(token)
        elif stack and stack_has_higher_precedence(stack, token):
            while stack and stack_has_higher_precedence(stack, token):
                queue.append(stack.pop())
            stack.append(token)
        else:
            stack.append(token)

    while stack:  # while stack is not empty
        queue.append(stack.pop())

    return queue


def eval_postfix(postfix) -> float:
    """
    Evaluate postfix expression to a single number
    """
    stack = []

    for token in postfix:
        if is_float(token):
            stack.append(float(token))
        elif token == 'i' or token == 'j':
            stack.append(complex(0.0, 1.0))
        elif token == 'pi':
            stack.append(math.pi)
        elif token == 'e':
            stack.append(cmath.exp(1))
        elif token == 'unary_minus':
            if len(stack) < 1:
                raise ValueError('Invalid expression')
            a = stack.pop()
            stack.append(-a)
        else:
            b = stack.pop()
            if token == '+':
                a = stack.pop()
                result = a + b
            elif token == '-':
                a = stack.pop()
                result = a - b
            elif token == '*':
                a = stack.pop()
                result = a * b
            elif token == '/':
                if b == 0:
                    raise ZeroDivisionError("Division by zero")
                a = stack.pop()
                result = a / b
            elif token == '^':
                a = stack.pop()
                result = a ** b
            elif token == 'sqrt':
                result = cmath.sqrt(b)
            elif token == 'sin':
                result = cmath.sin(b)
            elif token == 'cos':
                result = cmath.cos(b)
            elif token == 'tan':
                result = cmath.tan(b)

            stack.append(result)

    return stack.pop()


def evaluate(expression) -> float:
    return eval_postfix(infix_to_postfix(expression))


def truncate(number, decimals=0):
    if isinstance(number, complex):
        if decimals == 0:
            real_part = int(number.real)
            imag_part = int(number.imag)
        else:
            factor = 10.0 ** decimals
            real_part = int(number.real * factor) / factor
            imag_part = int(number.imag * factor) / factor

        # If the imaginary part is 0, return only the real part
        if imag_part == 0:
            return real_part
        else:
            return complex(real_part, imag_part)
    else:
        if decimals == 0:
            return int(number)
        factor = 10.0 ** decimals
        return int(number * factor) / factor


def apply_tolerance(value, tolerance=1e-8):
    """
    Fix approximation errors where a number is not recognized as 0 because it is e-18 or something.
    Works for both matrices and numbers
    """
    if isinstance(value, np.ndarray):
        return np.where(np.abs(value) < tolerance, 0, value)
    else:
        return 0 if np.abs(value) < tolerance else value


def lists_are_equal(l1, l2, tolerance=1e-8) -> bool:
    """
    Compare if two lists are equal while accounting for Python approximation errors
    """
    if len(l1) != len(l2):
        return False
    for val1, val2 in zip(l1, l2):
        if abs(val1 - val2) > tolerance:
            return False
    return True
