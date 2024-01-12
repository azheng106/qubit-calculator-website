from State import State


def get_constants_from_inputs(user_input):
    """
    Converts user input into a list of floats
    :param user_input: Input string from user
    :return: List of floats (input string separated by spaces)
    """
    try:
        return [float(x) for x in user_input.split()]
    except ValueError:
        return None


while True:
    input_str = input('Enter 8 constants c0 -> c7 separated by commas, or type "exit" to quit' + '\n')

    if input_str.lower() == 'exit' or input_str.lower() == 'quit':
        break

    constants = get_constants_from_inputs(input_str)

    if constants and len(constants) == 8:
        state = State(*constants)
        print("SLOCC Equivalence Class: " + str(state.get_state_class().name))  # .name removes the "EqClasses." prefix
    else:
        print("Invalid input. Please input 8 numeric values separated by spaces.")
