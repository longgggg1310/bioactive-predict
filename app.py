import numpy as np
from flask import Flask, request, render_template
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit import Chem
import pubchempy
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import wikipediaapi
from keras.models import load_model
import pickle

flask_app = Flask(__name__)
model2 = pickle.load(open("svm.pkl", "rb"))
model = load_model("network2.h5" )


@flask_app.errorhandler(404)
def not_found(error):
    return render_template('404.html'),404

@flask_app.route("/modal")
def modal():
    return render_template('modal.html', methods = ["POST"] )

@flask_app.route("/")
def Home():
    return render_template("index.html")

@flask_app.route("/predict", methods = ["POST"])
def predict():
    smile = ''.join([str(x) for x in request.form.values()])
    try:
        formula = getSmile(smile)
        name_text = getName(formula)
        mol = AllChem.MolFromSmiles(smile)
        AllChem.Compute2DCoords(mol)
        Draw.MolToFile(mol,"static/smiles.png",size=(400,300))
        full_name = change_name(smile)
        float_features = smiles_to_fp(smile)

        features1 = [np.array(float_features)]

        predictionML = model2.predict(features1)

        predictionML = returnType(predictionML)

        float_features= np.expand_dims(np.array(float_features), 0)

        print(float_features)
        predictionDL = model.predict(float_features)
        predictionDL= changeSmiles(predictionDL)
        # prediction= prediction.round(0).astype(int)
        print({'predictionDL': predictionDL, 'final_prediction': predictionDL, 'full_name': full_name, 'molecular': formula, 'predictionML': predictionML,
        'name_text': name_text})
        return render_template("result.html", predictionDL=predictionDL, predictionML = predictionML,
                               name=smile, full_name=full_name, molecular=formula, 
                               name_text=name_text)
    except Exception as e:
        print('Wrong fomular or wrong input')
        return not_found(404)
def getSmile(smile):
    mol = Chem.MolFromSmiles(smile)
    formula = CalcMolFormula(mol)
    return formula
def returnType(smile):
    if smile == 0 : 
        return 'ACTIVE'
    else: return 'INACTIVE'
def changeSmiles(smile) :
    if float(smile) >= 7:
        return 'ACTIVE'
    elif float(smile) <= 6:
        return 'INACTIVE'
    else: 
        return 'INTERMEDIATE'
def getName(smile):
    wiki_wiki = wikipediaapi.Wikipedia('en')
    page_py = wiki_wiki.page(smile)

    ("Page - Exists: %s" % page_py.exists())
    a = ("%s" % page_py.title)
    return a

def smiles_to_fp(smiles, method="maccs", n_bits=2048):
   

    # Convert smiles to RDKit mol object
    mol = Chem.MolFromSmiles(smiles)

    if method == "maccs":
        return np.array(MACCSkeys.GenMACCSKeys(mol))
    if method == "morgan2":
        return np.array(GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits))
    if method == "morgan3":
        return np.array(GetMorganFingerprintAsBitVect(mol, 3, nBits=n_bits))
    else:
        print(f"Warning: Wrong method specified: {method}." " Default will be used instead.")
        return np.array(MACCSkeys.GenMACCSKeys(mol))

def change_name(name):
    compounds = pubchempy.get_compounds(name, namespace='smiles')
    match = compounds[0]
    return match.iupac_name

if __name__ == "__main__":
    flask_app.run(debug=True)