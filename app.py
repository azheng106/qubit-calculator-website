from flask import Flask, request, jsonify, render_template

from util.evalutil import evaluate
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
        c0 = evaluate(data['c0'])
        c1 = evaluate(data['c1'])
        c2 = evaluate(data['c2'])
        c3 = evaluate(data['c3'])
        c4 = evaluate(data['c4'])
        c5 = evaluate(data['c5'])
        c6 = evaluate(data['c6'])
        c7 = evaluate(data['c7'])

        state = ThreeQubitState(c0, c1, c2, c3, c4, c5, c6, c7)

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
        if list(pairs)[-1] == '':
            pairs.popitem()  # Remove the last entry in dictionary if it is an empty pair ('': '')

        for key in pairs:  # Blank value entries default to 0
            if pairs[key] == '':
                pairs[key] = '0'
            if int(key) < 0:
                raise ValueError('Constant number < 0')

        qubit_value_pairs = {convert_to_basic_state(int(num_qubits), int(key)): evaluate(value)
                             for key, value in pairs.items()}

        # Fix approximation errors with Python sqrt
        tolerance = 1e-8
        keys = list(qubit_value_pairs.keys())
        values = list(qubit_value_pairs.values())
        for i in range(len(qubit_value_pairs) - 1):
            if abs(values[i] - values[i+1]) < tolerance:  # If values are close enough, then set equal to each other
                values[i+1] = values[i]
        qubit_value_pairs = dict(zip(keys, values))  # Recreate the dictionary with new values

        classification = is_entangled(qubit_value_pairs, int(num_qubits))

        return jsonify({'classification': classification.value})
    except Exception as e:
        return jsonify({'error': str(e)}), 400


if __name__ == '__main__':
    app.run(debug=False)
