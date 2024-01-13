from flask import Flask, request, jsonify, render_template

from State import State

app = Flask(__name__)


@app.route('/')
def index():
    """ Makes the website display when app.py is running and we visit http://localhost:5000"""
    return render_template('index.html')


@app.route('/calculate', methods=['POST'])
def calculate():
    """ Calculate and print the SLOCC class on the website when the user inputs 8 constants and presses 'submit' """
    data = request.json
    state = State(data['c0'], data['c1'], data['c2'], data['c3'],
                  data['c4'], data['c5'], data['c6'], data['c7'])
    classification = state.get_state_class().name
    return jsonify({'classification': classification})


if __name__ == '__main__':
    app.run(debug=True)
