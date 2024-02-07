from flask import Flask, request, jsonify, render_template
from util.evalutil import infix_to_postfix, eval_postfix
from util.state import State
from util.qubitutil import *

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/calculate/3qubit', methods=['POST'])
def calculate_slocc():
    try:
        """ Calculate and print the SLOCC class on the website when the user inputs 8 constants and presses 'submit' """
        data = request.json  # data (dictionary): {'c0': 'userinput', 'c1': 'userinput', ... , 'c7': 'userinput'}

        for key in data:
            if data[key] == '':  # Blank fields default to 0
                data[key] = '0'

        # Evaluate the expressions that the user inputted
        c0 = eval_postfix(infix_to_postfix(data['c0']))
        c1 = eval_postfix(infix_to_postfix(data['c1']))
        c2 = eval_postfix(infix_to_postfix(data['c2']))
        c3 = eval_postfix(infix_to_postfix(data['c3']))
        c4 = eval_postfix(infix_to_postfix(data['c4']))
        c5 = eval_postfix(infix_to_postfix(data['c5']))
        c6 = eval_postfix(infix_to_postfix(data['c6']))
        c7 = eval_postfix(infix_to_postfix(data['c7']))

        state = State(c0, c1, c2, c3, c4, c5, c6, c7)

        classification = state.get_state_class().value

        return jsonify({'classification': 'Belongs to ' + classification + ' SLOCC class'})
    except Exception as e:
        return jsonify({'error': str(e)}), 400  # Error code 400


@app.route('/calculate/nQubit', methods=['POST'])
def calculate_entanglement():
    try:
        data = request.json  # Ex. of data: {'n': '4', 'pairs': {'3': '123', '2': '13', '4': '12', '': ''}}
        num_qubits = data['n']

        pairs = data['pairs']
        pairs.popitem()  # Remove the last entry in dictionary, which is an empty pair ('': '')
        for key in pairs:  # Blank value entries default to 0
            if pairs[key] == '':
                pairs[key] = '0'
        qubit_value_pairs = {convert_to_basic_state(int(num_qubits), int(key)): eval_postfix(infix_to_postfix(value)) for key, value in pairs.items()}

        #print(f"Original state: {qubit_value_pairs}")
        classification = 'State is unknown'
        if is_entangled(qubit_value_pairs, int(num_qubits)):
            classification = 'State is entangled'

        return jsonify({'classification': classification})
    except Exception as e:
        return jsonify({'error': str(e)}), 400


if __name__ == '__main__':
    app.run(debug=False)
