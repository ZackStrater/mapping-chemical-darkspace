
import sys
sys.path.append('../')
from src.find_fragments import find_fragment, fragmentize
from src.fragments_library import special_cases, biomolecules, peptide_amino_acids, heterocycles, \
    common_aromatic_heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments

libraries = [heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments, special_cases]


from flask import Flask, render_template, request, jsonify
app = Flask(__name__)

@app.route('/', methods=['GET'])
def index():
    return render_template('molecule_analyzer.html')


@app.route('/get_fragments', methods=['POST'])
def get_fragments():
    user_data = request.json
    input_string = user_data['input_string']
    fragments, molecule = fragmentize(input_string, *libraries)
    return jsonify({'fragments': fragments})


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, debug=True)



