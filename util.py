import math
import re
import cmath
from enum import Enum


class Precedence(Enum):
    # Lower values indicate higher precedence
    FUNCTION = 0
    POWER = 1
    MULTIPLY = 2
    ADD = 3
    OPEN_PARENTHESIS = 99  # When comparing any operator to '(', '(' is always lower precedence


precedences = {
    '(': Precedence.OPEN_PARENTHESIS.value,
    'sin': Precedence.FUNCTION.value,
    'cos': Precedence.FUNCTION.value,
    'tan': Precedence.FUNCTION.value,
    'sqrt': Precedence.FUNCTION.value,
    '^': Precedence.POWER.value,
    '*': Precedence.MULTIPLY.value,
    '/': Precedence.MULTIPLY.value,
    '+': Precedence.ADD.value,
    '-': Precedence.ADD.value,
}


def stack_has_higher_precedence(stack, current_token) -> bool:
    """
    Return True if the operator on top of the stack has >= precedence than the current token, False otherwise
    """
    stack_precedence = precedences.get(stack[-1], float('inf'))
    current_precedence = precedences.get(current_token, float('inf'))
    return stack_precedence <= current_precedence


def is_negative_sign(tokens, index):
    """
    Determine if the "-" is a negative number (True) or minus sign (False)
    :param index: index of the "-"
    :param tokens: list of tokens
    """
    if index == 0 or tokens[index - 1] in '*/+-^(':
        return True
    return False


def infix_to_postfix(infix) -> list:
    """
    Use Shunting Yard Algorithm to translate mathematical infix to postfix
    https://www.youtube.com/watch?v=Wz85Hiwi5MY
    :param infix: Input mathematical expression
    :return: Postfix
    """
    tokens = re.findall(r'[\d.]+|[+\-*/()^]|sqrt|i|j|sin|cos|tan', infix)  # infix = '1+1' => tokens = ['1', '+', '1']
    stack = []  # Operator stack
    queue = []  # Output queue

    for index, token in enumerate(tokens):
        if token.isnumeric() or token in ['i', 'j']:  # Numbers are automatically pushed into stack
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
        if token.isnumeric():
            stack.append(float(token))
        elif token == 'i' or token == 'j':
            stack.append(complex(0.0, 1.0))
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
