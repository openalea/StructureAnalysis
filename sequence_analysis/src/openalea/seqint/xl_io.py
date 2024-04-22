"""
.. module:: xl_io
   :platform: Unix, Windows
   :synopsis: Ecriture de fichiers xl.

"""
from mmap import mmap, ACCESS_READ
from xlrd import open_workbook
from xlutils.copy import copy


def XLWriteStates(inputfile, outputfile, seq, rvar, skipnwfreq=False):
    """Ecrit la sequence d'etats dans un fichier excel.

    Args:
        inputfile: nom du fichier XL d'origine
        outputfile: nom du fichier XL de sortie
        seq: sequences d'etats et d'observations
        rvar: variable correspondant a readmode dans seq
        skipnwfreq: si True, ignore les lignes de wfreq < 0

   .. note::
        variable 0 correspond a la variable etat

   """
    book = open_workbook(inputfile, formatting_info=True)
    assert(book.nsheets == 1)

    # rs = rb.sheet_by_index(0)
    cpbook = copy(book)

    # feuille a copier
    rs = book.sheets()[0]

    # feuille a ecrire
    ws = cpbook.get_sheet(0)

    # insere deux colonnes
    col = ws.col(rs.ncols)
    col = ws.col(rs.ncols+1)

    # colonne readmode
    readmode = [i for i in range(rs.ncols) if "READMODE" == str(rs.cell(0,i).value)][0]

    # colonne wfreq
    wfreq = [i for i in range(rs.ncols) if "WFREQ" == str(rs.cell(0,i).value)][0]

    # colonne isfirst
    fst = [i for i in range(rs.ncols) if "ISFIRST" == str(rs.cell(0,i).value)][0]

    # colonne islast
    lst = [i for i in range(rs.ncols) if "ISLAST" == str(rs.cell(0,i).value)][0]

    ws.write(0, rs.ncols, "STATES")

    ws.write(0, rs.ncols+1, "READMODE_CTRL")

    var = rvar
    last_state = 0 # dernier etat lu
    s = 0 # sequence
    t = 0 # position dans la sequence
    j = 2 # ligne

    for i in range(1, rs.nrows):
        if t == 0:
            # trouve la longueur de la sequence, wreq < 0 exclus ou pas
            pred_py_length = 0 # longueur theorique de la sequence python
            xl_length = 0 # longueur de la sequence xl
            j = i+xl_length # ligne courante
            while (j < rs.nrows) and len(str(rs.cell(j,lst).value)) > 0 \
                and (int(rs.cell(j,lst).value) == 0):
                if ((float(rs.cell(j,wfreq).value) > -1) or (not(skipnwfreq))):
                    pred_py_length += 1
                j += 1
            fin_fichier = len(str(rs.cell(j,lst).value)) == 0
            if (j < rs.nrows) and not(fin_fichier) \
                and ((float(rs.cell(j,wfreq).value) > -1) or (not(skipnwfreq))): # cellule islast
                pred_py_length += 1
            if (j < rs.nrows) and not(fin_fichier):
                j += 1
            xl_length = j - i
            # augmenter i en consequence
            if (pred_py_length < 3) and not(fin_fichier) and (len(seq[s]) < 3):
                # on a oublie de supprimer la sequence python ?
                s += 1
        if (pred_py_length > 2):
            if ((float(rs.cell(i,wfreq).value) <= -1) and skipnwfreq):
                # saute la valeur et reecrit l'etat precedent
                ws.write(i, rs.ncols, last_state)
                # ecrit readmode pour controle
                ws.write(i, rs.ncols+1, int(rs.cell(i,readmode).value))
            else:
                state = seq[s][t][0]
                assert(int(rs.cell(i,readmode).value) == seq[s][t][var])
                ws.write(i, rs.ncols, state)
                # ecrit readmode pour controle
                ws.write(i, rs.ncols+1, seq[s][t][var])
                last_state = state
                t += 1
            if(int(rs.cell(i,lst).value) == 1):
                # derniere valeur lue
                s += 1
                t = 0
        else:
            # recopie du fichier excel avec etat -1
            while (i <= j) and (i < rs.nrows):
                if len(str(rs.cell(i,readmode).value)) > 0:
                    ws.write(i, rs.ncols, -1)
                    # ecrit readmode pour controle
                    ws.write(i, rs.ncols+1, int(rs.cell(i,readmode).value))
                i += 1

    cpbook.save(outputfile)

def XLFindTransition(inputfile, edeb, efin, var=None):
    """Trouve toutes les sequences qui contiennent la transition
       de edeb vers efin pour la variable var (par defaut, STATES)
    """

    book = open_workbook(inputfile, formatting_info=True)
    assert(book.nsheets == 1)
    i = 0

    # feuille a lire
    rs = book.sheets()[0]

    # colonne subj
    scol = [i for i in range(rs.ncols) if "SUBJ" == str(rs.cell(0,i).value)][0]

    # colonne text_no
    text_no = [i for i in range(rs.ncols) if "TEXT_NO" == str(rs.cell(0,i).value)][0]

    # colonne isfirst
    fst = [i for i in range(rs.ncols) if "ISFIRST" == str(rs.cell(0,i).value)][0]

    # colonne islast
    lst = [i for i in range(rs.ncols) if "ISLAST" == str(rs.cell(0,i).value)][0]

    # colonne states
    states = [i for i in range(rs.ncols) if "STATES" == str(rs.cell(0,i).value)][0]

    iv = var
    if iv is None:
        iv = states

    # recherche de transition
    i = 1
    tr = [] # textes contenant la transition

    etat_courant = -1
    etat_precedent = -1
    fin_fichier = False

    while (i < rs.nrows) and not(fin_fichier):
        etat_courant = int(rs.cell(i,iv).value)
        if int(rs.cell(i,fst).value) == 1:
            # nouveau texte
            etat_precedent = -1
        if etat_courant > -1:
            if (etat_courant == efin) and \
                (etat_precedent == edeb):
                    tr += [(int(rs.cell(i,scol).value), \
                        int(rs.cell(i,text_no).value))]
            etat_precedent = etat_courant
        i += 1
        fin_fichier = len(str(rs.cell(i,lst).value)) == 0

    tr = set(tr)
    return tr
