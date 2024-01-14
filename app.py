from flask import Flask, request, jsonify, render_template
from util import infix_to_postfix, eval_postfix
from State import State

app = Flask(__name__)


@app.route('/')
def index():
    """ Makes the website display when app.py is running and we visit http://localhost:5000"""
    return render_template('index.html')


@app.route('/calculate', methods=['POST'])
def calculate():
    """ Calculate and print the SLOCC class on the website when the user inputs 8 constants and presses 'submit' """
    data = request.json  # data (dictionary): {'c0': 'userinput', 'c1': 'userinput', ... , 'c7': 'userinput'}

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

    classification = state.get_state_class().name
    return jsonify({'classification': classification})


if __name__ == '__main__':
    app.run(debug=True)
