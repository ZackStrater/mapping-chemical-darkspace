
import sys
sys.path.append('../')
from src.find_fragments import find_fragment, fragmentize
from src.fragments_library import special_cases, biomolecules, peptide_amino_acids, heterocycles, \
    common_aromatic_heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments

libraries = [heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments, special_cases]
names, molecule = fragmentize(r'IC1=CNC(=O)C=C1', *libraries)
print(names)




from flask import Flask, render_template, request, jsonify
app = Flask(__name__)

@app.route('/', methods=['GET'])
def index():
    return render_template('molecule_analyzer.html')


@app.route('/get_fragments', methods=['POST'])
def get_fragments():
    user_data = request.json
    a = user_data['a']
    root_1 = _solve_quadratic(a)
    return jsonify({'root_1': root_1})


def _solve_quadratic(a):
    return 'molecule' + a


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, debug=True)



