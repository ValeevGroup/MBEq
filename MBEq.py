from pdb import set_trace
import copy

print "\n\n\n"
print "                                            Welcome to Automatic Equation Generator\n\n"
print "                                                  Liguo Kong and Marcel Nooijen\n"
print "                               Department of Chemistry, University of Waterloo, Waterloo, ON, Canada\n\n\n"




# notation:
#         'sf': 'spin-free'
#         'so': 'spin-orbital'
#         'km': 'MR normal order by Mukherjee and Kutzelnigg'






#**************
#  Directory  *
#**************




Public    =   ['Classes and Related',
               'Index and Excitation Operator Base',
               'Contraction Functions',
               'KM-GWT',
               'Tools Library',
               'Spatial Orbital Canonicalization: connected',
               'Spatial Orbital Canonicalization: permutation',
               'Spin Orbital Canonicalization',
               'spin expansion',
               'SRCC',
               'MR-F12',
               'Factorization',
               'Latex',
               'Ohio  Tensor Contraction Engine Interface']








#******************************
#  System Control Parameters  *
#******************************






system_switch = {             'use spin orbital' : 0,                 # this global variable controls whether we are doing spin or spatial orbital work, which will affect
                                                                      # , e.g., coeff(), or even generalIndex, we also need completely different merge_all/canonical functions
                                                                      # at the beginning, this is due to the demand of Ohio TCE people, but we might use it in our own work later.
                                                                      # WHENEVER DOING SPIN ORBITAL WORK, TURN THIS ON
                                  'normal order' : 'T',
                                                                      # if system_switch['normal order'] == 'T',   use traditional normal order contraction scheme
                                                                      # if system_switch['normal order'] == 'KM',  use the KM`-normal order contraction scheme; spin-free version shall not be used
                                                                      # if system_switch['normal order'] == 'V',   genuine vacuum normal order, based on tweaking of KM normal order
                                                                      # MN normal ordering is deprecated
                             'matrix value table': 'r12',
                                 'explicit spin' : 0,          
                              'canonical scheme' : 'p',               # how to canonicalize:  'p' means 'permu', works generally; 'c' means 'connect', works for spin-free connected expression not involving 
                                                                      # quantities of rank higher than 2
                         'GWT CU truncate rank' : 1,                  # In MR normal order expansion, truncate cumulant rank; if N = 0, no truncation, keep all;
                   'EXCLUDE Particle Ind from CU': 1,                 # if reference WFN is a CAS reference function, then cumulants containing particle indices will be zero;
                                                                      # if it is a general wavefunction, then cumulants containing particle indices may not be zero                                                                                                       
                     'Allow Internal Contraction': 0,                 # either 0 or 1; used for MR normal order expansion; this is useful for expanding DM into cumulants, via a trick                                                                      
                                 'expand tensor' : [],                # this is the contraction variable used in MN/MN2-normal order theory, e.g. ['DEN', 'Cbar']                                       
                                    'Max DM rank': 2,
                                 'DM coefficient': 'normal'} 
                                                                      # in this definition of density matrix, D_ij^kl is not <L|k^+l^+ ji |R>, so use this global variable to
                                                                      # control; if D_ij^kl = <L|k^+l^+ ji |R>, we set system_switch['DM coefficient'] == 'normal'
                                                                      # otherwise it is '1/n!'
                                                                      # it should be 'normal 'when we do KM theory ! because that is how we define the lists








#**************************************
#       Parameters used in the code   *
#**************************************






typelist = ['p', 'e', 'h', 'm', 'g']


type_value = [0, 1, 2, 3, 4]


dm_or_cumu_names = ['Lambda', 'lambda', 'lambda1', 'Lambda1', 'DEN', 'Gamma', 'gamma', 'd', 'd1', 'd2', 'd3', 'd4']

dic_type_value = {'P': 0,
                  'p': 1,              # this capital letter stands for alpha spin orbital index; used in explict spin cases       
                  'E': 10,
                  'e': 11,
                  'H': 100,
                  'h': 101,
                  'M': 1000,
                  'm': 1001,
                  'g': 10001,
                  'A': 0,
                  'a': 1,
                  'I': 100,
                  'i': 101}



dic_active_value = {-1: 3,
                     0: 31,
                     1: 301}



dic_fix_value = {0: 7,
                 1: 700}



dic_matname_value_all = {'f': 41,
                    'nf': 41,
                   'kmf': 41,
                     'v': 43,
                    'nv': 43,
                   'kmv': 43,
                  'lvsA': 47,
                  'lvsB': 53,
                  'lvsC': 59,
                  'lvsD': 61,
                     't': 79,
                    't1': 83,
                    't2': 97,
                   'ltA': 101,
                    'lc': 103,
                     's': 107,
                    's1': 109,
                    's2': 113,
                    's3': 119,
                   'lsB': 121,
                    's4': 127,
                   'lsA': 131,
                    's5': 133,
                    's6': 134,
                     'w': 137,
                   'lsC': 139,
                    'nd': 199,
                     'X': 142,         # 'X' will not appear in tenta_final()
                     'z': 100000,      # 'z' is used occasionally in mcscf
                    'lx': 150,         # 'lx' is used occasionally in mrucc  
                 'delta': 151,
                  'F_t1': 155,
                  'F_t2': 157,
                'F(t1)_': 159,
                'F(t2)_': 161,
                 'hbar1': 169,         # for ccsd_h
                 'hbar2': 179,
                    'g1': 171,         # for steom_g
                    'g2': 173,
                    'g3': 177,
                     'M': 183,
                   'DEN': 187,
                  'Cbar': 193,
                    't3': 199,
                   'eta': 209,
                 'gamma': 217,
                 'Gamma': 211,
                'lambda': 221,
                'Lambda': 212,
               'lambda1': 223,
               'Lambda1': 227,
                     'd': 233,
                    'd1': 239,
                    'd2': 240,
                    'd3': 241,
                    'd4': 242,
                     'G': 243,
                     'K': 244,
                 'D_INV': 245,
              'Cbar_INV': 246,
                    'x1': 253,         # this quantity is used to MN2-n.o.; I will also use this to temperarily simulate lambda1_inv to avoid introducing new mats
                    'x2': 257,
                    'x3': 259,
                     'r': 263}         # this corresponds to the r12 integral


dic_matname_value_r12 = {'t2': 1,
                  'Gamma': 10,
                 'Lambda': 20,
                      'f': 30,
                      'v': 40,
                'Lambda1': 50,
                      'r': 60,
                     'x1': 70,      # this corresponds to the r12 integral
                     'd1': 80,
                     'd2': 90,
                     'd3': 100,
                  'delta': 110}


dic_matname_value_srsd = {'f': 41,
                     'v': 43,
                    't1': 83,
                    't2': 97,
                 'delta': 151,
                    't3': 199}



# this is to allow more freedom; it only affects the looking of the final expressions.
if system_switch['matrix value table'] == 'r12': dic_matname_value = dic_matname_value_r12
elif system_switch['matrix value table'] == 'srsd': dic_matname_value = dic_matname_value_srsd
elif system_switch['matrix value table'] == 'showcase': dic_matname_value = dic_matname_value_all 
else: 
    print "please specify"
    set_trace()


new_dic_mat_value = {'t2': 1,
                      'd': 5}



dic_matname_va = { 't': 't',
                  't1': 't',
                  't2': 't',
                  't3': 't',
                  's1': 's',
                  's2': 's',
                  's3': 's',
                  's4': 's',
                  's5': 's',
                  's6': 's',
                   'w': 'w',
                   'f': 'f',
                  'nf': 'F',
                 'kmf': 'F',
                   'v': 'v',
                  'nv': 'V',
                 'kmv': 'V',
                   'd': 'd',
                  'd1': 'd',
                  'd2': 'd',
                  'd3': 'd',
                  'd4': 'd',
                  'nd': 'nd',
                   'M': 'M',
                   't': 't',
                'lvsA': '\wp',
                'lvsB': '\\aleph',
                'lvsC': '\Im',
                'lvsD': '\Re',
                'ltA': '\mho',
                  'lc': '\sigma',
                 'lsA': '\zeta',
                 'lsB': '\\xi',
                 'lsC': '\eta',
                   'X': 'X',
                  'lx': '\Delta',
               'delta': '\delta',
                   'I': 'I',
                 'cor': ' E ',
                'hbar1': '\\bar{H}',
                'hbar2': '\\bar{H}',
                   'g1': 'g',
                   'g2': 'g',
                   'g3': 'g',
                   'g4': 'g',
                   'J1': 'J_eom_1',
                   'J2': 'J_eom_2',
                    'A': 'A',
                    'B': 'B',
                    'C': 'C',
                    'D': 'D',
                    'E': 'E',
                    'F': 'F',
                  'DEN': '\\gamma',
                   'x1': '\\xi',
                   'x2': '\\xi',
                   'x3': '\\xi',
                    'G': 'G',
                    'K': 'K',
                'D_INV': '{\\bar{\mathbb{D}}}',
             'Cbar_INV': '{\\bar{\\eta}}',
                 'Cbar': '\\bar{D}',
                  'eta': '\\eta',
                'Gamma': '\\Gamma',
               'Lambda': '\\Lambda',
              #'Lambda1': '\\Lambda',
              'Lambda1': '\\Gamma',
               'lambda': '\\lambda',
              'lambda1': '\\lambda',
                'gamma': '\\gamma',
                    #'r': '\\bar{r}'}
                    'r': 'R'}















#--------------------------------------------------------------------------------------------------------------------------------------

#                                          Section 1:         CLASSES AND RELATED METHODS

#--------------------------------------------------------------------------------------------------------------------------------------









#**************************************************
#                class-Index                      *
#**************************************************





class Index:
    # in the following line, num = "" can be simplified to "", but the form below is more reminicent
    def __init__(self, type = "", num = "", att = "", fix = 0, active = 0, other = ["",    "",     ""]):
        self.type = type
        self.num = num
        self.att = att
        self.fix = fix
        self.active = active
        self.other = other

    def equal(one, two):
        if one.type == two.type and one.num == two.num and one.fix == two.fix and one.active == two.active: return 1
        else: return 0


    def similar(one, two):
        if one.type == two.type: return 1
        else: return 0


    def inArray(ind, array):
        for a in array:
            if a.equal(ind): return 1
        return 0


    def gvIndexNo(ind, array):
        AA = range(len(array))
        for a in AA:
            if array[a].equal(ind): return a
        return res  # -1000 indicates: not in array


    def gvType(indx):
        return indx.type


    def gvIndex(index):
        return index.type + index.num


    def printindex(index):
        print index.gvIndex()






def null_index_att_s(terms):           # refer to the above
    result = []
    for xx in terms:
        result.append(null_index_att(xx))
    return result








#**************************************************
#              class-Uoperator                    *
#**************************************************





class Uoperator:
    def __init__(self, upperIndicees = [], lowerIndicees = [], cha = -1, other = []):
        self.upperIndicees = upperIndicees
        self.lowerIndicees = lowerIndicees
        self.cha = cha
        self.other = other


    def gvUppInd(self):
        return self.upperIndicees


    def gvLowInd(self):
        return self.lowerIndicees


    def equals(self, operator2):
        result = 0
        set1 = makeIndexSet(self.upperIndicees, self.lowerIndicees)
        set2 = makeIndexSet(operator2.upperIndicees, operator2.lowerIndicees)
        for a in set1.indexPairs:
            if(a.inArray(set2.indexPairs)): result += 1
            else: pass
        if(result == len(set1.indexPairs) == len(set2.indexPairs)): return 1
        else: return 0









#**************************************************
#              class-TypeContract                 *
#**************************************************





class TypeContract:

    def __init__(self, indexList1 = [], indexList2 = []):

        """ indexList1 has the index no. of labels of first operator and indexList2 has of second operator, so they record which indices are contracted """
        self.indexList1 = indexList1
        self.indexList2 = indexList2





#**************************************************
#               class-Contract                    *
#**************************************************





class Contract:

    """ class of sets of contractions has several particle-type contractions and several hole-type contractions """
    def __init__(self, pTypeContract = TypeContract(), hTypeContract = TypeContract()):
        self.pTypeContract = pTypeContract
        self.hTypeContract = hTypeContract





#**************************************************
#            class-ContableArrays                 *
#**************************************************





class ContableArrays:

    """ the class of arrays of labels which can be contracted. e.g.: a hole label can not be contracted in a p-type contraction """

    def __init__(self, array11 = [], array12 = [], array21 = [], array22 = []):
        self.array11 = array11
        self.array12 = array12
        self.array21 = array21
        self.array22 = array22




#**************************************************
#              class-MatElement                   *
#**************************************************





class MatElement:

    """ class for matrix elements, mainly attributes: name, matLowerIndicees, matUpperIndicees """
    def __init__(self, name = "", matUpperIndicees = [], matLowerIndicees = [], connectivity = 0, other = []):
        self.name = name
        self.matLowerIndicees = matLowerIndicees
        self.matUpperIndicees = matUpperIndicees
        self.connectivity = connectivity                               # used in STORNGLY CONNECTED COUPLED CLUSTER. I think this
        self.other = other                                             # attribute is only used temporarily, to select strongly connected terms


    def equals(element1, element2):
        """ MatElement class method: determine whether element1 and element2 are the same or not """ 
        result = 0
        if(element1.name == element2.name):
            if(len(element1.matUpperIndicees) != len(element2.matUpperIndicees)):
                result = 0
            else:
                set1 = makeIndexSet(element1.matUpperIndicees, element1.matLowerIndicees)
                set2 = makeIndexSet(element2.matUpperIndicees, element2.matLowerIndicees)
                for a in set1.indexPairs:
                    if(a.inArray(set2.indexPairs)):
                        result = 1
                    else:
                        result = 0
                        break
        return result


    def mat_2_indexset(self):

        """ convert a 'mat' object to an IndexSet object; if later on we do conversion back, then mat.connectivity and mat.other information is lost """
        indp = []
        for n in range(len(self.matLowerIndicees)):
            indp += [IndexPair(self.matUpperIndicees[n], self.matLowerIndicees[n])]
        return IndexSet(indp)  



def find_outer(mats):
    """ find external indices """
    if len(mats) == 0: return Uoperator([], [])
    elif len(mats) == 1: return Uoperator(mats[0].matUpperIndicees, mats[0].matLowerIndicees)
    else:
        yy = mats[0]
        for xx in range(len(mats)-1):
            yy = find_outer_2([yy, mats[xx + 1]])
        return Uoperator(yy.matUpperIndicees, yy.matLowerIndicees)


def find_outer_2(mat2):
    """ the difference from find_outer_2_bk is that here by the minor change we avoid the problem happening with, say, diagnaol mat element f^a_a, not tested yet """
    upp = copy.deepcopy(mat2[0].matUpperIndicees + mat2[0].matLowerIndicees)
    low = copy.deepcopy(mat2[1].matUpperIndicees + mat2[1].matLowerIndicees)
    L1 = []
    L2 = []
    LL = []
    for xxx in range(len(upp)):
        if(not upp[xxx].inArray(low)): L1 += [upp[xxx]]
        else: low = delete(low, upp[xxx].gvIndexNo(low))
    LL = L1 + low
    n = len(LL)                        # a simple check here
    if(n - 2*(n/2) != 0):
        print "\n\n Odd number of overall indices in find_outer_2()!"
        set_trace
    if(len(LL) > 0):
        x = len(LL)/2
        L1 = LL[:x]
        L2 = LL[x:]
    return MatElement("", L1, L2, mat2[0].connectivity, mat2[0].other)











#**************************************************
#              class-Coefficient                  *
#**************************************************





class Coefficient:

    """ two attributes: const, matElement, e.g. for the term 0.5 t^ab_ij E^ab_ij, Coefficient = 0.5 t^ab_ij, Coefficient.const = 0.5, Coefficient.matElement = [t^ab_ij] """
    def __init__(self, const = 1, matElement = []):
        self.const = const
        self.matElement= matElement


    def equals(c1, c2):
        result = 0
        co1 = []
        co2 = []
        if(len(c2.matElement) == len(c1.matElement)):
            for a in range(len(c1.matElement)):
                for b in range(len(c2.matElement)):
                    if(c1.matElement[a].equals(c2.matElement[b])):
                        if(b in co1): pass
                        else:
                            co1.append(b) 
                            result += 1
            if(result == len(c1.matElement)): return 1
            else: return 0
        return 0


    def construct_0_op_term(self):
        return Term(coefficient = self)



def coeff_s_times_coeff_s(cos1, cos2):
    return [Coefficient(const = a.const * b.const, matElement = a.matElement + b.matElement) for a in cos1 for b in cos2]


def coeff_s_to_2_op_terms(cos):
    return [Term(coefficient = ss) for ss in cos]
    






#**************************************************
#                   class-Term                    *
#**************************************************





class Term:
    def __init__(self, coefficient = Coefficient(), uOperators = [Uoperator([], [])], type = -1, other = []):     # add 'other' for future generalization
        self.coefficient = coefficient                                            
        self.uOperators = uOperators
        self.type = type
        self.other = other

    def gvFirstOp(term):
        return term.uOperators[0]

    def gvSecondOp(term):
        return term.uOperators[1]

    def equals(t1, t2):
        result = 1
        if(t1.coefficient.equals(t2.coefficient)): pass
        else: result = 0
        for a in range(len(t1.uOperators)):
            if(t1.uOperators[a].equals(t2.uOperators[a])): pass
            else: result = 0
        return result

    def inArray(term, terms):
        res = 0
        for a in terms:
            if(a.equals(term)): 
                res = 1
                break
        return res


    def inArray2(term, terms):
        res = 0
        pos = 0
        for nn, a in enumerate(terms):
            if(a.equals(term)): 
                res = 1
                pos = nn 
        return [res, pos]

    def gvIndexNo(term, terms):
        res = 0
        for a in range(len(terms)):
            if(terms[a].equals(term)):
                res = a
                break
            else:
                pass
        return res

    def null_term_type(term):
        term.type =-1
        return term


    def reorder_term(term1):
        term = copy.deepcopy(term1)
        term = reorder(term)
        for i in range(len(term.coefficient.matElement)):
            term = term.reorder_mat_i(i)
        return term


    def null_index_att(term):
        """ reset index.att to default """
        for inde in term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees:
            inde.att = ""
        for xx in term.coefficient.matElement:
            for yy in xx.matUpperIndicees + xx.matLowerIndicees:
                yy.att = ""
        return term



    def null_index_fix(term):              # reset index.fix to default
        for inde in term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees:
            inde.fix = 0
        for xx in term.coefficient.matElement:
            for yy in xx.matUpperIndicees + xx.matLowerIndicees:
                yy.fix = 0
        return term


  
    def null_index_active(term):           # reset index.active to default
        for inde in term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees:
            inde.active = 0
        for xx in term.coefficient.matElement:
            for yy in xx.matUpperIndicees + xx.matLowerIndicees:
                yy.active = 0
        return term


    def reorder_mat_i(term1, i):                                             
        """ arrange indexpairs in term1.coefficient.matElement[i] in an convenient order, pp-pe-...-gg """
        mattotal = []
        term = copy.deepcopy(term1)
        for k in range(25):
            mattotal = mattotal + [[]]
        map = []
        totalnum = []
        typelist = ["p", "e", "h", "m", "g"]
        for dd in typelist:
            for ff in typelist:
                map = map + [[dd, ff]]
        oo = term.coefficient.matElement[i]
        for ii in range(len(oo.matUpperIndicees)):
            jj = [oo.matLowerIndicees[ii].type, oo.matLowerIndicees[ii].type]
            mattotal[map.index(jj)] = mattotal[map.index(jj)] + [ii]
        upp = []
        low = []
        for kk in mattotal:
            for ll in kk:
                upp = upp + [oo.matUpperIndicees[ll]]
                low = low + [oo.matLowerIndicees[ll]]
        oo.matUpperIndicees = upp
        oo.matLowerIndicees = low
        term.coefficient.matElement[i] = MatElement(term.coefficient.matElement[i].name, upp, low, term.coefficient.matElement[i].connectivity, term.coefficient.matElement[i].other)
        return term


    def reorder_mats(self, preferorder):
        """ preferorder gives the ones you want to appear first """
        mats1 = []
        mats2 = []
        if len(preferorder) == 1:
            for xx in self.coefficient.matElement:
                if xx.name == preferorder[0]: mats1.append(xx)
                else: mats2.append(xx)
            return Term(Coefficient(self.coefficient.const, mats1 + mats2), self.uOperators)
        else: set_trace()  # only 1 is implemented so far
       

    def hasmat_name(self, matname):
        """ tell whether termm contain a mat named 'matnmae'; return [yesorno (0/1), numbers] """
        find = 0
        num = 0
        for xx in self.coefficient.matElement: 
            if xx.name == matname:
                find = 1
                num += 1
        return [find, num]



    def hasmat(self, matt):
        """ tell whether termm contain matt or not; if yes, return 'pos'; otherwise return -1 """
        res = -1
        nummats = len(self.coefficient.matElement)
        for xx in range(nummats): 
            if(matt.equals(self.coefficient.matElement[xx])):
                res = xx
                break
        return res

    def hasmats(self, mats):
        """ tell if term has all mats; and separate them into 2 groups """
        termm = copy.deepcopy(self)
        hasall = 1
        poss = []
        mm = copy.deepcopy(self.coefficient.matElement)
        for xx in mats:
            pos = self.hasmat(xx)
            if pos == -1:
               hasall = 0
               break
            else: poss.append(pos)
        if hasall:
            poss.sort()
            poss.reverse()
            mm = copy.deepcopy(self.coefficient.matElement)
            for xx in poss: mm.pop(xx)      # remove mats from the term
        return [hasall, mm]


    def has_indices_more_than_twice(self):
        find = 0
        allinds = give_all_ind(self)
        uniqueinds = give_indlist(self)
        for xx in uniqueinds:
            if appear_time(xx, allinds) > 2:
                xx.printindex()
                find = 1
                break
        return find


    def give_ex_deex_rank(self):
        """ Give the excitation and deexcitation rank in diagrammatic sense.
    
            assume no 'general' indices exist """
        ex = 0
        deex = 0
        for yy in self.uOperators[0].upperIndicees:
            if yy.type in ['p', 'e', 'P', 'E']:
                ex += 1
            elif yy.type in ['h', 'm', 'H', 'M']:
                deex += 1
            else:
                print "general indices? Not handled. Break it first"
                printterm(term)
                set_trace()
        for yy in self.uOperators[0].lowerIndicees:
            if yy.type in ['p', 'e', 'P', 'E']:
                deex += 1
            elif yy.type in ['h', 'm', 'H', 'M']:
                ex += 1
            else:
                print "general indices? Not handled. Break it first"
                printterm(term)
                set_trace()
        return [ex, deex, ex-deex]



    def replace_2b_cumu_with_dm_1(self):
        if(system_switch['use spin orbital'] == 1): CumulantName = 'lambda'
        else: CumulantName = 'Lambda'
        flag = 0
        others = []
        flag2 = 1
        for xx in self.coefficient.matElement:
            if(xx.name == CumulantName and len(xx.matUpperIndicees) == 2 and flag2 == 1):     # since we only replace ONE cumulant, we use flag2 to monitor
                flag += 1
                eta = xx
                flag2 = 0
            else:
                others += [xx]
        if(flag == 0):
            printterm(self)
            pdb.set_trace()
            return
        elif(system_switch['use spin orbital'] == 1):    
            return [Term(Coefficient(self.coefficient.const, others + [MatElement('gamma', eta.matUpperIndicees, eta.matLowerIndicees)]), self.uOperators),
                    Term(Coefficient(-1.0*self.coefficient.const, others + [MatElement('lambda1', eta.matUpperIndicees[:1], eta.matLowerIndicees[:1]), MatElement('lambda1', eta.matUpperIndicees[1:], eta.matLowerIndicees[1:])]), self.uOperators),
                    Term(Coefficient(self.coefficient.const, others + [MatElement('lambda1', eta.matUpperIndicees[:1], eta.matLowerIndicees[1:]), MatElement('lambda1', eta.matUpperIndicees[1:], eta.matLowerIndicees[:1])]), self.uOperators)]
        else:
            return [Term(Coefficient(self.coefficient.const, others + [MatElement('Gamma', eta.matUpperIndicees, eta.matLowerIndicees)]), self.uOperators),
                    Term(Coefficient(-1.0*self.coefficient.const, others + [MatElement('Lambda1', eta.matUpperIndicees[:1], eta.matLowerIndicees[:1]), MatElement('Lambda1', eta.matUpperIndicees[1:], eta.matLowerIndicees[1:])]), self.uOperators),
                    Term(Coefficient(self.coefficient.const * 0.5, others + [MatElement('Lambda1', eta.matUpperIndicees[:1], eta.matLowerIndicees[1:]), MatElement('Lambda1', eta.matUpperIndicees[1:], eta.matLowerIndicees[:1])]), self.uOperators)]



    def replace_2b_dm_with_cumu_1(self):
        term = self
        if system_switch['use spin orbital'] == 1: DMName = 'gamma'
        else: DMName = 'Gamma'
        flag = 0
        others = []
        flag2 = 1
        for xx in term.coefficient.matElement:
            if  flag2 == 1 and xx.name == DMName and len(xx.matUpperIndicees) == 2:     # since we only replace ONE DM, we use flag2 to monitor
                flag += 1
                eta = xx
                flag2 = 0
            else: others += [xx]
        if flag == 0:
            printterm(term)
            set_trace()
        elif system_switch['use spin orbital'] == 1:    
            return [Term(Coefficient(term.coefficient.const, others + [MatElement('lambda', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators),
                    Term(Coefficient(term.coefficient.const, others + [MatElement('lambda1', eta.matUpperIndicees[:1], eta.matLowerIndicees[:1]), MatElement('lambda1', eta.matUpperIndicees[1:], eta.matLowerIndicees[1:])]), term.uOperators),
                    Term(Coefficient(-1.0 * term.coefficient.const, others + [MatElement('lambda1', eta.matUpperIndicees[:1], eta.matLowerIndicees[1:]), MatElement('lambda1', eta.matUpperIndicees[1:], eta.matLowerIndicees[:1])]), term.uOperators)]
        else:
            return [Term(Coefficient(term.coefficient.const, others + [MatElement('Lambda', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators),
                    Term(Coefficient(term.coefficient.const, others + [MatElement('Lambda1', eta.matUpperIndicees[:1], eta.matLowerIndicees[:1]), MatElement('Lambda1', eta.matUpperIndicees[1:], eta.matLowerIndicees[1:])]), term.uOperators),
                    Term(Coefficient(term.coefficient.const * (-0.5), others + [MatElement('Lambda1', eta.matUpperIndicees[:1], eta.matLowerIndicees[1:]), MatElement('Lambda1', eta.matUpperIndicees[1:], eta.matLowerIndicees[:1])]), term.uOperators)]




    def replace_12b_sf_dm_with_delta_1(self):
        if(system_switch['use spin orbital'] != 0): set_trace()
        term = self
        DMNames = ['Gamma', 'Lambda1']
        flag = 0
        others = []
        flag2 = 0
        for xx in term.coefficient.matElement:
            if flag == 0 and xx.name in DMNames and len(xx.matUpperIndicees) in [1, 2]:     # since we only replace ONE DM, we use flag2 to monitor
                flag = 1
                ga = xx
            else:
                others += [xx]
        if(flag == 0):
            printterm(term)
            set_trace()
            return
        else:
            if len(ga.matUpperIndicees) == 2:
                return [Term(Coefficient(term.coefficient.const * 4.0, others + [MatElement('delta', ga.matUpperIndicees[:1], ga.matLowerIndicees[:1]),MatElement('delta', ga.matUpperIndicees[1:], ga.matLowerIndicees[1:])]), term.uOperators),
                        Term(Coefficient(term.coefficient.const * (-2.0), others + [MatElement('delta', ga.matUpperIndicees[:1], ga.matLowerIndicees[1:]),MatElement('delta', ga.matUpperIndicees[1:], ga.matLowerIndicees[:1])]), term.uOperators)]
            elif len(ga.matUpperIndicees) == 1:
                return [Term(Coefficient(term.coefficient.const * 2.0, others + [MatElement('delta', ga.matUpperIndicees, ga.matLowerIndicees)]), term.uOperators)]
            else: set_trace()



    def replace_sf_3_dm_with_cumu_1(term):
        """ Replace spin-free 3-RDM by 1/2-cumulant (spin-free), as in the MolPhys paper by Kutzelnigg, Eq. (75)  """
        flag = 0
        others = []
        flag2 = 1
        for xx in term.coefficient.matElement:
            if(xx.name == 'd3' and flag2 == 1):  
                flag += 1
                gamma3 = xx
                flag2 = 0
            else: others += [xx]
        if(flag == 0):
            printterm(term)
            set_trace()
            return
        elif(system_switch['use spin orbital'] == 0):    
            [P1, P2, P3, Q1, Q2, Q3] = gamma3.matUpperIndicees +  gamma3.matLowerIndicees 
            Coeffs = [  Coefficient(1.0,  [MatElement('Lambda1', [P1], [Q1]),MatElement('Lambda1', [P2], [Q2]), MatElement('Lambda1', [P3], [Q3]) ]),
                        Coefficient(0.25, [MatElement('Lambda1', [P1], [Q2]),MatElement('Lambda1', [P2], [Q3]), MatElement('Lambda1', [P3], [Q1]) ]),
                        Coefficient(0.25, [MatElement('Lambda1', [P1], [Q3]),MatElement('Lambda1', [P2], [Q1]), MatElement('Lambda1', [P3], [Q2]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P1], [Q2]),MatElement('Lambda1', [P2], [Q1]), MatElement('Lambda1', [P3], [Q3]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P1], [Q1]),MatElement('Lambda1', [P2], [Q3]), MatElement('Lambda1', [P3], [Q2]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P1], [Q3]),MatElement('Lambda1', [P2], [Q2]), MatElement('Lambda1', [P3], [Q1]) ]),
                        Coefficient(1.0,  [MatElement('Lambda1', [P1], [Q1]),MatElement('Lambda', [P2, P3], [Q2, Q3]) ]),
                        Coefficient(1.0,  [MatElement('Lambda1', [P2], [Q2]),MatElement('Lambda', [P1, P3], [Q1, Q3]) ]),
                        Coefficient(1.0,  [MatElement('Lambda1', [P3], [Q3]),MatElement('Lambda', [P1, P2], [Q1, Q2]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P1], [Q2]),MatElement('Lambda', [P2, P3], [Q1, Q3]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P1], [Q3]),MatElement('Lambda', [P2, P3], [Q2, Q1]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P2], [Q1]),MatElement('Lambda', [P1, P3], [Q2, Q3]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P2], [Q3]),MatElement('Lambda', [P1, P3], [Q1, Q2]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P3], [Q1]),MatElement('Lambda', [P1, P2], [Q3, Q2]) ]),
                        Coefficient(-0.5, [MatElement('Lambda1', [P3], [Q2]),MatElement('Lambda', [P1, P2], [Q1, Q3]) ])]
            TempTerm = Term(Coefficient(term.coefficient.const, others), term.uOperators) 
            res = []
            for xx in Coeffs: res.append(Term(Coefficient(TempTerm.coefficient.const * xx.const, others + xx.matElement), TempTerm.uOperators))
            return res











def select_ex_deex_rank(terms, rank):
    """ rank is a list of three numbers [ex rank, de rank, ex - de]; -1 means 'not compare this number' """
    res = []
    for yy in terms:
        yyrank =  yy.give_ex_deex_rank()
        okay = 1
        if (yyrank[0] == 2 * rank[0] or rank[0] == -1) and  (yyrank[1] == 2 * rank[1] or rank[1] == -1) and (yyrank[2] == 2 * rank[2] or rank[2] == -1): 
            res.append(yy)
    return res




def terms_has_indices_more_than_twice(terms):
    """ if this happens, something is wrong """
    rr = range(len(terms))
    for xx in rr:
        if terms[xx].has_indices_more_than_twice(): 
            print xx
            printterm(terms[xx])
            set_trace()
    return 
        

def terms_reorder_mats(terms, preferorder):
    res = []
    for xx in terms: res.append(xx.reorder_mats(preferorder))
    return res
       


def null_terms_type(terms):
    result = []
    for xx in terms:
        result += [xx.null_term_type()]
    return result

def reorder_term_s(terms):
    res = []
    for xx in terms:
        res += [xx.reorder_term()]
    return res

def null_index_fix_s(terms):           # refer to the above
    result = []
    for xx in terms:
        result += [null_index_fix(xx)]
    return result




def null_index_active_s(terms):        # refer to the above
    result = []
    for xx in terms:
        result += [null_index_active(xx)]
    return result


             


#**************************************************
#                  class-Loop                     *
#**************************************************





class Loop:

        """ a loop has a starting point, an ending point and the middle points """

	def __init__(self, type = "", startP = 0, endP = 0, middles = []):
		self.type = type
		self.startP = startP
		self.endP = endP
		self.middles = middles





#**************************************************
#              class-IndexPair                    *
#**************************************************





class IndexPair:

    def __init__(self, upperIndex, lowerIndex, ordering = 1):
        self.upperIndex = upperIndex
        self.lowerIndex = lowerIndex
        self.ordering = ordering

	

    def show(pair):
        return " (" + pair.upperIndex.gvIndex() + ")" + "(" + pair.lowerIndex.gvIndex() + ") "



    def similar(p1, p2):
        """ tells whether the two pairs of labels are similar or not """
        if(p1.upperIndex.type == p2.upperIndex.type and p1.lowerIndex.type == p2.lowerIndex.type):
            return 1
        else:
            return 0

	
    def gvIndexNo(indP, arr):
  
        """ gives the position of the index pair "indP" in the array "arr" """
        res = 0
        for a in range(len(arr)):
            if(arr[a].upperIndex.equal(indP.upperIndex) and arr[a].lowerIndex.equal(indP.lowerIndex)):
                res = a
                break
            else:
                pass
        return res



    def rmEl(pair, array):

        """ removes the index pair object "pair" from the array of index pair objects "array" """
        return array[:pair.gvIndexNo(array)] + array[pair.gvIndexNo(array) + 1:]



    def inArray(indP, arr):

        """ tells whether the indexpair object "indP" is present in the array "arr" or not """
        res = 0
        for a in arr:
            if(a.upperIndex.equal(indP.upperIndex) and a.lowerIndex.equal(indP.lowerIndex)):
                res = 1
                break
            else:
                pass
        return res





#**************************************************
#              class-IndexSet                     *
#**************************************************





class IndexSet:

    """ contains all the index pairs of a single matrix element or an operator """
    def __init__(self, indexPairs = []):
        self.indexPairs = indexPairs
        
     
    def get_list(self):
        A = []
        B = []
        for n in range(len(self.indexPairs)):
            A += [self.indexPairs[n].upperIndex]
            B += [self.indexPairs[n].lowerIndex]
        return A + B



"""  a few IndexSet-related functions """



def show(self):
    Indicees([self]).show()



def findfixind(ind):

    """ find fixed indices from 'ind' """

    result =   [ [],  [],  [],  [],  [],  [],  [],  [],  []]
    typelist = ["p", "e", "h", "m", "g", "P", "E", "H", "M"]
    for ab in range(len(ind.indexSets)):
        for bc in range(len(ind.indexSets[ab].indexPairs)):
            if(ind.indexSets[ab].indexPairs[bc].upperIndex.fix == 1):
                if(ind.indexSets[ab].indexPairs[bc].upperIndex.inArray(result[typelist.index(ind.indexSets[ab].indexPairs[bc].upperIndex.type)])== 0):
                    result[typelist.index(ind.indexSets[ab].indexPairs[bc].upperIndex.type)] += [ind.indexSets[ab].indexPairs[bc].upperIndex]
            if(ind.indexSets[ab].indexPairs[bc].lowerIndex.fix == 1):
                if(ind.indexSets[ab].indexPairs[bc].lowerIndex.inArray(result[typelist.index(ind.indexSets[ab].indexPairs[bc].lowerIndex.type)])== 0):
                    result[typelist.index(ind.indexSets[ab].indexPairs[bc].lowerIndex.type)] += [ind.indexSets[ab].indexPairs[bc].lowerIndex]

    return result





#**************************************************
#              class-Indicees                     *
#**************************************************





class Indicees:

      """ this class is important for index manipulation; practically it is a list of 'IndexSet' """
      def __init__(self, indexSets = []):
            self.indexSets = indexSets



      def renumber(ind):

            """ renumber indices in 'Ind' which is an object of Indicees. fixed indices will not be reordered, so first we pick out fixed indices, group them """
            """ and exclude their numbers from being used in following renumbering procedure; otherwise confusion raised """

            res = Indicees()
            res.indexSets = []
            diction = []
            dKeys = []
            PP = []
            EE = []
            HH = []
            MM = []
            GG = []
            P2 = []
            E2 = []
            H2 = []
            M2 = []
            G = []
            Fix = findfixind(ind)
            for x in Fix[0] + Fix[5]:
                PP.append(int(x.num))
            for x in Fix[1] + Fix[6]:
                EE.append(int(x.num))
            for x in Fix[2] + Fix[7]:
                HH.append(int(x.num))
            for x in Fix[3] + Fix[8]:
                MM.append(int(x.num))
            for x in Fix[4]:
                GG.append(int(x.num))
            P_max = 0
            E_max = 0
            H_max = 0
            M_max = 0
            G_max = 0
            if PP:                    # make sure the list is not empty
                P_max = max(PP)
            if EE:
                E_max = max(EE)
            if HH:
                H_max = max(HH)
            if MM:
                M_max = max(MM)
            if GG:
                G_max = max(GG)       # the above 5 nines are newly added; these variables are used to make sure that after renumbering, 
                                      # the 'num's of all non-fixed indices are higher than any of the 'num's of the fixed indices for each type
            for j in range(1000):
                i = j + 1
                if(i not in PP and i > P_max):
                    P2 += [i]
                if(len(P2) == 40):
                    break
            for j in range(1000):
                i = j + 1
                if(i not in EE and i > E_max):
                    E2 += [i]
                if(len(E2) == 40):
                    break
            for j in range(1000):
                i = j + 1
                if(i not in HH and i > H_max):
                    H2 += [i]
                if(len(H2) == 40):
                    break
            for j in range(1000):
                i = j + 1
                if(i not in MM and i > M_max):
                    M2 += [i]
                if(len(M2) == 40):
                    break
            for j in range(1000):
                i = j + 1
                if(i not in GG and i > G_max):
                    G += [i]
                if(len(G) == 20):
                    break
            p2 = 0
            e2 = 0
            h2 = 0
            m2 = 0
            g = 0
            p2Count = P2[p2]
            e2Count = E2[e2]
            h2Count = H2[h2]
            m2Count = M2[m2]
            gCount = G[g]
            for a in range(len(ind.indexSets)):
                  set = IndexSet()
                  set.indexPairs = []
                  for b in range(len(ind.indexSets[a].indexPairs)):
                        yyx = ind.indexSets[a].indexPairs[b].upperIndex
                        if(yyx.fix == 1):
                              set.indexPairs += [IndexPair(yyx, Index(), 1)]
                        elif(yyx.inArray(dKeys)):
                              set.indexPairs += [IndexPair(diction[yyx.gvIndexNo(dKeys)], Index(), 1)]
                        else:
                              dKeys += [yyx]
                              if(yyx.type in ["P", "p"]): 
                                    set.indexPairs += [IndexPair(Index(yyx.type, str(p2Count), yyx.att, 0, yyx.active, yyx.other), Index(), 1)]  
                                    diction += [Index(yyx.type, str(p2Count), yyx.att, 0, yyx.active, yyx.other)]
                                    p2 += 1                                               # there is no conflict. 
                                    p2Count = P2[p2]
                              elif(dKeys[-1].type in ["H", "h"]):
                                    set.indexPairs += [IndexPair(Index(yyx.type, str(h2Count), yyx.att, 0, yyx.active,  yyx.other), Index(), 1)]
                                    diction += [Index(dKeys[-1].type, str(h2Count),  yyx.att, 0,  yyx.active,  yyx.other)]
                                    h2 += 1
                                    h2Count = H2[h2]
                              elif(dKeys[-1].type in ["M", "m"]):
                                    set.indexPairs += [IndexPair(Index(yyx.type, str(m2Count),  yyx.att, 0,  yyx.active,  yyx.other), Index(), 1)]
                                    diction += [Index(yyx.type, str(m2Count), yyx.att, 0, yyx.active, yyx.other)]
                                    m2 += 1
                                    m2Count = M2[m2]
                              elif(dKeys[-1].type in ["E", "e"]):
                                    set.indexPairs += [IndexPair(Index(yyx.type, str(e2Count), yyx.att, 0, yyx.active, yyx.other), Index(), 1)]
                                    diction += [Index(yyx.type, str(e2Count), yyx.att, 0, yyx.active, yyx.other)]
                                    e2 += 1
                                    e2Count = E2[e2]
                              elif(dKeys[-1].type == "g"):
                                    set.indexPairs += [IndexPair(Index("g", str(gCount), yyx.att, 0, yyx.active, yyx.other), Index(), 1)]
                                    diction += [Index(yyx.type, str(gCount), yyx.att, 0, yyx.active, yyx.other)]
                                    g += 1
                                    gCount = G[g]
                              else:
                                    print "\ntype unexpected\n"
                                    yyx.printindex()
                                    set_trace()
                  for b in range(len(ind.indexSets[a].indexPairs)):
                        yyx = ind.indexSets[a].indexPairs[b].lowerIndex
                        if(yyx.fix == 1):
                              set.indexPairs[b].lowerIndex = yyx
                        elif(yyx.inArray(dKeys)):
                              set.indexPairs[b].lowerIndex = diction[yyx.gvIndexNo(dKeys)]
                        else:
                              dKeys += [yyx]
                              if(yyx.type in ["P", "p"]):
                                    set.indexPairs[b].lowerIndex = Index(yyx.type, str(p2Count), yyx.att, 0, yyx.active, yyx.other)
                                    diction += [Index(yyx.type, str(p2Count), yyx.att, 0, yyx.active, yyx.other)]
                                    p2 += 1
                                    p2Count = P2[p2]
                              elif(yyx.type in ["H", "h"]):
                                    set.indexPairs[b].lowerIndex = Index(yyx.type, str(h2Count), yyx.att, 0, yyx.active, yyx.other)
                                    diction += [Index(yyx.type, str(h2Count), yyx.att, 0, yyx.active, yyx.other)]
                                    h2 += 1
                                    h2Count = H2[h2]
                              elif(yyx.type in ["M", "m"]):
                                    set.indexPairs[b].lowerIndex = Index(yyx.type, str(m2Count), yyx.att, 0, yyx.active, yyx.other)
                                    diction += [Index(yyx.type, str(m2Count), yyx.att, 0, yyx.active, yyx.other)]
                                    m2 += 1
                                    m2Count = M2[m2]
                              elif(yyx.type in ["E", "e"]):
                                    set.indexPairs[b].lowerIndex = Index(yyx.type, str(e2Count), yyx.att, 0, yyx.active, yyx.other)
                                    diction += [Index(yyx.type, str(e2Count), yyx.att, 0, yyx.active, yyx.other)]
                                    e2 += 1
                                    e2Count = E2[e2]
                              elif(yyx.type == "g"):
                                    set.indexPairs[b].lowerIndex = Index("g", str(gCount), yyx.att, 0, yyx.active, yyx.other)
                                    diction += [Index(yyx.type, str(gCount), yyx.att, 0, yyx.active, yyx.other)]
                                    g += 1
                                    gCount = G[g]
                              else:
                                    print "\ntype unexpected\n"
                                    yyx.printindex()
                                    set_trace()
                  res.indexSets += [copy.deepcopy(set)]
            return res



      def show(self):
          str = ""
          for a in range(len(self.indexSets)):
	        str = str + "["
	        for b in self.indexSets[a].indexPairs:
		      str = str + " ( " + b.upperIndex.gvIndex() + " " + b.lowerIndex.gvIndex() + " )"
	        str = str + " ] "
          print str







#**************************************************
#               Constructions                     *
#**************************************************





def index_construction(indexstr):

    """ construct an index; indexstr = e.g. ['p', '1', {}]; {} = e.g. {'fix': 0, 'active': 0}, which is for fine tuning; if no fine tuning, the third argument {} must be supplied """

    if(len(indexstr) != 3):
        print "\n error in function index_construction@eq_lib \n"
        set_trace()
    indtype = indexstr[0]
    indnum  = indexstr[1]
    IXX = Index()
    fixx = IXX.fix
    activee = IXX.active
    attt = IXX.active
    if('att' in indexstr[2].keys()):
        attt = indexstr[2]['att']
    if('fix' in indexstr[2].keys()):
        fixx = indexstr[2]['fix']
    if('active' in indexstr[2].keys()):
        activee = indexstr[2]['active']
    return Index(indtype, indnum, attt, fixx, activee) 



def indices_construction(indsstr):

    """ construct a list of indices; insstr is a list of 'indexstr', e.g., ['p', '1', {}, 'h', '2', {}], as in the previous function """

    result = []
    if(len(indsstr) == 0):
        pass
    elif(len(indsstr)%3 != 0):         # '3' arguments are needed to construct one index, as in the previous function
        print "\n error in indices_contruction@eq_lib \n"
        set_trace()
    else:
        index_quantity = len(indsstr)/3
        for xx in range(index_quantity):
            #pdb.set_trace()
            result += [index_construction(indsstr[3 * xx: 3 * xx + 3])]
    return result




def op_construction(upp_low_indlist_str):

    """ construct one matrixelement from [indsstr, indsstr]  """

    if(len(upp_low_indlist_str) != 2):
        print "\n error in op_construction \n"
        set_trace()
    uppinds = indices_construction(upp_low_indlist_str[0])
    lowinds = indices_construction(upp_low_indlist_str[1])
    return Uoperator(uppinds, lowinds)



def ops_construction(ops_str):

    """ construct a list of operators from ops_str; ops_str = e.g. [upp_low_indlist_str, upp_low_indlist_str]  """

    result = []
    for xx in ops_str:
        result += [op_construction(xx)]
    return result



def one_mat_construction(upp_low_indlist_str):

    """ construct one matrixelement from [name, indsstr, indsstr]  """

    if(len(upp_low_indlist_str) != 3):
        print "\n error in one_mat_construction \n"
        set_trace()
    namee = upp_low_indlist_str[0]
    uppinds = indices_construction(upp_low_indlist_str[1])
    lowinds = indices_construction(upp_low_indlist_str[2])
    return MatElement(namee, uppinds, lowinds)



def mats_construction(matsindsstr):
 
    """ construct a list of matrixelements; matsindsstr = e.g., [upp_low_indlist_str, upp_low_indlist_str] """

    result = []
    for xxx in matsindsstr:
        result += [one_mat_construction(xxx)]
    return result



def coefficient_construction(co_str):

    """ construct one Coefficient object; co_str = e.g., [2.0, matsindsstr]; if there is no matrixelement, please supply '[]' to matsindsstr  """

    return Coefficient(co_str[0], mats_construction(co_str[1]))



def term_construction(term_str):
   
    """ construct a term; term_str = e.g. [2.0, matsindsstr, ops_str] """

    return Term(coefficient_construction(term_str[:2]), ops_construction(term_str[2]))



def term_with_one_mat_construction(term_str):

    """ construct a term which has one matrixelement and whose indices in matlement and operator correspond, e.g., the frequently used H/F/V/T terms; term_str = [const, name, uppstr, dowstr], the same as in 
        function one_mat_construction """
 
    mats = [one_mat_construction(term_str[1:])]                                # the first argument is the 'const'   
    opp = Uoperator(mats[0].matUpperIndicees, mats[0].matLowerIndicees)
    return Term(Coefficient(term_str[0], mats), [opp])


    


def direct_op_construction(upp_low_indlist):

    if(len(upp_low_indlist)%2 != 0):
        print "\n error \n"
        set_trace()
    if(len(upp_low_indlist) == 0):
        return Uoperator([], [])
    else:
        return Uoperator(upp_low_indlist[:len(upp_low_indlist)/2], upp_low_indlist[len(upp_low_indlist)/2:])



def direct_ops_construction(list_of_inds):
    """ Construct a list of operators from list_of_inds, which is e.g. [ [p1, h1], [p2, h2]]. """

    result = []
    if len(list_of_inds) == 0: return [Uoperator([], [])]
    else:
        for xx in list_of_inds: result += [direct_op_construction(xx)]
        return result



def direct_one_mat_construction(name_inds):
    """ Construct one matrixelement from [name, inds].  """

    if(len(name_inds) != 2):
        print "\n error in one_mat_construction \n"
        pdb.set_trace()
    return MatElement(name_inds[0], name_inds[1][:len(name_inds[1])/2], name_inds[1][len(name_inds[1])/2:])



def direct_mats_construction(list_of_name_inds):
    """ Construct a list of matrixelements; matsindsstr = e.g., [ [name, inds], [name, inds] ]. """

    result = []
    for xxx in list_of_name_inds:
        result += [direct_one_mat_construction(xxx)]
    return result



def direct_coefficient_construction(co_str):

    """ construct one Coefficient object; co_str = e.g., [2.0, list_of_name_inds]; if there is no matrixelement, please supply '[]' to the second argument  """

    return Coefficient(co_str[0], direct_mats_construction(co_str[1]))



def direct_term_construction(term_str):
   
    """ construct a term; term_str = e.g. [2.0, matsindsstr, ops_str] """

    return Term(direct_coefficient_construction(term_str[:2]),    \
                    direct_ops_construction(term_str[2]))





    



def showindlist(list):                 # print index list
    res = ""
    for ind in list: res += ind.gvIndex() + ", "
    print res



def showOperator(operator):
    if len(operator.upperIndicees) == 0: return ""
    else:
        str = "E ( "
        for a in operator.upperIndicees:
            str += a.gvIndex() + " "
        str += ") ( "
        for b in operator.lowerIndicees:
            str = str + b.gvIndex() + " "
        str += ") "
        return str


def show_num_operator(operator):
    if(len(operator.upperIndicees) == 0):
        strr = ""
    else:
        strr = "E( "
        for a in operator.upperIndicees:
            strr = strr + str(a)
        strr = strr + ";"
        for b in operator.lowerIndicees:
            strr = strr + str(b)
        strr = strr + ")"
    return strr


def show_num_operator_s(operators):
    strr = ""
    for xx in operators:
        strr += show_num_operator(xx) + " "
    return strr 

def print_num_operator_s(operators):
    print show_num_operator_s(operators)


def printOperator(operator):
    print showOperator(operator)

def print_num_operator(operator):
    print show_num_operator(operator)

def print_num_operator_s(operators):
    print show_num_operator_s(operators) 

def print_num_operator_s_s(operators_s):
    for xx in operators_s:
        print_num_operator_s(xx)
    return



def showMatElement(matElement):
    """ returns the matrix element in string form """
    str = matElement.name + " ( "
    for a in matElement.matUpperIndicees:
        str = str + a.gvIndex() + " "
    str = str + ") ( "
    for b in matElement.matLowerIndicees:
        str = str + b.gvIndex() + " "
    str = str + ") "
    return str



def printMatElements(mats):
    S = ""
    for xx in mats: S += showMatElement(xx)
    print S


def printMatElement(matElement):
    print showMatElement(matElement)





def showTerm(te):
    """ returns the term in string form """
    string = "[ " + str(te.coefficient.const) + " ] "
    for b in te.uOperators: string += showOperator(b)
    for a in te.coefficient.matElement: string = string + showMatElement(a)
    return string





def showterms(terms):
    """ returns the terms (array of objects of class term) in string form"""
    res = "   "
    for a in terms: 
        if terms.index(a) < len(terms) - 1: res += showTerm(a) + "\n" + " + "
        else: res += showTerm(a)
    return res



def printterms(terms):
    print showterms(terms)



def printterm(term):
    printterms([term])










#--------------------------------------------------------------------------------------------------------------------------------------

#                                              Index and Excitation Operator Base

#--------------------------------------------------------------------------------------------------------------------------------------






class Index_OP_Base_Separator:
    def __init__(self, start = "Y"):
        self.start = type







#*****************************
#        Create Index        *
#*****************************



def pfactory(n = 500):
    i = 1
    res = []
    while(i<n+1):
        x = Index("p", str(i))
        res.append(x) 
        i += 1
    return res





def hfactory(n = 500):
    i = 1
    res = []
    while(i<n+1):
        x = Index("h", str(i))
        res.append(x)
        i += 1
    return res



def efactory(n = 500):
    i = 1
    res = []
    while(i<n+1):
        x = Index("e", str(i), "", 0, 1)                                     # '1' here denotes 'active'
        res.append(x)
        i += 1
    return res



def mfactory(n = 20):
    i = 1
    res = []
    while(i<n+1):
        x = Index("m", str(i), "", 0, 1)                                     # '1' here denotes 'active'
        res.append(x)
        i += 1
    return res



def gfactory(n = 20):
    i = 1
    res = []
    while(i<n+1):
        x = Index("g", str(i))
        res.append(x)
        i += 1
    return res





    
#*******************************************
#       Spin-orbital Excitation Operator   *
#*******************************************                                   
        
        
    
def so_EX_d():                              # spin orbital operator; so the coefficient becomes 1/n!
    res = []
    i = 0
    HF = hfactory()
    PF = pfactory()
    if 1:
        h1 = HF()[i]
        h2 = HF[i + 1]
        p1 = PF[i]
        p2 = PF[i + 1]
        dmat = MatElement("t2", [h1, h2], [p1, p2])
        dco = Coefficient(0.25, [dmat])
        dop = Uoperator([p1, p2], [h1, h2])
    return [Term(dco, [dop])]
    


def so_EX_sd():
    res = []
    i = 0
    HF = hfactory()
    PF = pfactory()
    if 1:
        h1 = HF[i]
        h2 = HF[i+1]
        p1 = PF[i]
        p2 = PF[i + 1]
        smat = MatElement("t1", [h1], [p1])
        sco = Coefficient(1.0, [smat])
        sop = Uoperator([p1], [h1])
        dmat = MatElement("t2", [h1, h2], [p1, p2])
        dco = Coefficient(0.25, [dmat])
        dop = Uoperator([p1, p2], [h1, h2])
    return [Term(sco, [sop]), Term(dco, [dop])]


def so_EX_t():                         # give ONE triple excitation operator
    HF = hfactory()
    PF = pfactory()
    res = []
    i = 0
    h1 = HF[i]
    h2 = HF[i + 1]
    h3 = HF[i + 2]
    p1 = PF[i]
    p2 = PF[i + 1]
    p3 = PF[i + 2]
    mat = MatElement("t3", [h1, h2, h3], [p1, p2, p3])
    co = Coefficient(1.0, [mat])  # 1/36 coefficient will be taken care in so_triple_coeff_s()     
    op = Uoperator([p1, p2, p3], [h1, h2, h3])
    return [Term(co, [op])]



        
#***********************************************************
#      Hamiltoniann and Transformed Hamiltonian Operator   *
#***********************************************************




def so_H():
    opa = Uoperator([gfactory()[0]], [gfactory()[1]])
    opb = Uoperator([gfactory()[0],gfactory()[1]], [gfactory()[2], gfactory()[3]])
    mat1 = MatElement("f", [gfactory()[1]],[gfactory()[0]])
    mat2 = MatElement("v", [gfactory()[2], gfactory()[3]], [gfactory()[0], gfactory()[1]])
    co1 = Coefficient(1.0, [mat1])
    co2 = Coefficient(0.25, [mat2])
    fterm = Term(co1, [opa])
    vterm = Term(co2, [opb])
    return [fterm] + [vterm]




def so_H_particular():                 # this is a special spin orbital H for speical purpose, it differs from so_H by the coefficient
    opa = Uoperator([gfactory()[0]], [gfactory()[1]])
    opb = Uoperator([gfactory()[0],gfactory()[1]], [gfactory()[2], gfactory()[3]])
    mat1 = MatElement("f", [gfactory()[1]],[gfactory()[0]])
    mat2 = MatElement("v", [gfactory()[2], gfactory()[3]], [gfactory()[0], gfactory()[1]])
    co1 = Coefficient(1.0, [mat1])
    co2 = Coefficient(0.5, [mat2])
    fterm = Term(co1, [opa])
    vterm = Term(co2, [opb])
    return [fterm] + [vterm]






def H():                               # return Hamiltonian in the form of list of terms
    opa = Uoperator([gfactory()[0]], [gfactory()[1]])
    opb = Uoperator([gfactory()[0],gfactory()[1]], [gfactory()[2], gfactory()[3]])
    mat1 = MatElement("f", [gfactory()[1]],[gfactory()[0]])
    mat2 = MatElement("v", [gfactory()[2], gfactory()[3]], [gfactory()[0], gfactory()[1]])
    co1 = Coefficient(1.0, [mat1])
    co2 = Coefficient(0.5, [mat2])
    fterm = Term(co1, [opa])
    vterm = Term(co2, [opb])
    return [fterm] + [vterm]


   

def ccsd_h():                          # Here we define \bar{H} = e^{-T} H e^T as intermediates, assuming that only
                                       # up to 2-body operators are used
                                       # Probably rasing the 'num' of indices often needed.
    opa = Uoperator([gfactory()[0]], [gfactory()[1]])
    opb = Uoperator([gfactory()[0],gfactory()[1]], [gfactory()[2], gfactory()[3]])
    mat1 = MatElement("hbar1", [gfactory()[1]],[gfactory()[0]])
    mat2 = MatElement("hbar2", [gfactory()[2], gfactory()[3]], [gfactory()[0], gfactory()[1]])
    co1 = Coefficient(1.0, [mat1])
    co2 = Coefficient(0.5, [mat2])
    fterm = Term(co1, [opa])
    vterm = Term(co2, [opb])
    return [fterm] + [vterm]





def sf_km_H():                         # return Hamiltonian in the form of list of terms
    opa = Uoperator([gfactory()[0]], [gfactory()[1]])
    opb = Uoperator([gfactory()[0],gfactory()[1]], [gfactory()[2], gfactory()[3]])
    mat1 = MatElement("kmf", [gfactory()[1]],[gfactory()[0]])
    mat2 = MatElement("kmv", [gfactory()[2], gfactory()[3]], [gfactory()[0], gfactory()[1]])
    co1 = Coefficient(1.0, [mat1])
    co2 = Coefficient(0.5, [mat2])
    fterm = Term(co1, [opa])
    vterm = Term(co2, [opb])
    return [fterm] + [vterm]





def normal_H():                        
    opa = Uoperator([gfactory()[0]], [gfactory()[1]])
    opb = Uoperator([gfactory()[0],gfactory()[1]], [gfactory()[2], gfactory()[3]])
    mat1 = MatElement("nf", [gfactory()[1]],[gfactory()[0]])
    mat2 = MatElement("nv", [gfactory()[2], gfactory()[3]], [gfactory()[0], gfactory()[1]])
    co1 = Coefficient(1.0, [mat1])
    co2 = Coefficient(0.5, [mat2])
    fterm = Term(co1, [opa])
    vterm = Term(co2, [opb])
    return [fterm] + [vterm]


















#--------------------------------------------------------------------------------------------------------------------------------------

#                                                Contraction Functions

#--------------------------------------------------------------------------------------------------------------------------------------










class ContractionSeparator:
    def __init__(self, start = "Y"):
        self.start = type










def makeContableArrays(fullTerm):

        """ obtain the arrays that can be contracted;  this method needs to be changed when more types of labels are introduced """

	ary11 = fullTerm.gvFirstOp().gvUppInd()
	ary12 = fullTerm.gvSecondOp().gvUppInd()
	ary21 = fullTerm.gvFirstOp().gvLowInd()
	ary22 = fullTerm.gvSecondOp().gvLowInd()
	array11 = []
	array12 = []
	array21 = []
	array22 = []
	for a in range(len(ary11)):
		if(ary11[a].gvType() == "p" or ary11[a].gvType() == "e"):
			pass
		else:
			array11 = array11 + [a]
		if(ary21[a].gvType() == "h" or ary21[a].gvType() == "m"):
			pass
		else:
			array21 = array21 + [a]
	for b in range(len(ary12)):
		if(ary12[b].gvType() == "h" or ary12[b].gvType() == "m"):
			pass
		else:
			array12 = array12 + [b]
		if(ary22[b].gvType() == "p" or ary22[b].gvType() == "e"):
			pass
		else:
			array22 = array22 + [b]
	return ContableArrays(array11, array12, array21, array22)





def tContractions(k, list1, list2):

        """ makes all possible k-contractions from two arrays; these array are just the positions of labels in the operators """

	array1 = combi(k, list1)
	array2 = permu(k, list2)
	contractions = []
	for a in array1:
		for b in array2:
			type = TypeContract(a,b)
			contractions = contractions + [type]
	return contractions



def contractions(k, arrays):

        """ gives all possible k-contractions from a ContableArrays object; gives an array of objects of the class Contract """

	n = min(len(arrays.array11), len(arrays.array22))
	p = min(len(arrays.array21), len(arrays.array12))
	contracts = []
	if(k <= n and k <= p):
		for a in range(k+1):
			pContrac = tContractions(k-a, arrays.array21, arrays.array12)
			nContrac = tContractions(a, arrays.array11, arrays.array22)
			if(len(pContrac) == 0):
				pContrac = [TypeContract()]
			if(len(nContrac) == 0):
				nContrac = [TypeContract()]
			for b in range(len(pContrac)):
				for c in range(len(nContrac)):
					contraction = Contract(pContrac[b], nContrac[c])
					contracts = contracts + [contraction]
		return contracts
	elif(k > min(p,n) and k <= max(p,n)):
		for a in range(min(p,n)+1):
			if(p>n):			
				pContractions = tContractions(k-a, arrays.array21, arrays.array12)
				nContractions = tContractions(a, arrays.array11, arrays.array22)
			else:
				pContractions = tContractions(a, arrays.array21, arrays.array12)
				nContractions = tContractions(k-a, arrays.array11, arrays.array22)
			if(len(pContractions) == 0):
				pContractions = [TypeContract()]
			if(len(nContractions) == 0):
				nContractions = [TypeContract()]
			for b in pContractions:
				for c in nContractions:
					contraction = Contract(b, c)
					contracts = contracts + [contraction]
		return contracts
	elif(k > max(p,n) and k <= p+n):
		for a in range(p+n-k+1):
			if(p>n):			
				pContractions = tContractions(p-a, arrays.array21, arrays.array12)
				nContractions = tContractions(k-p+a, arrays.array11, arrays.array22)
			else:
				pContractions = tContractions(k-n+a, arrays.array21, arrays.array12)
				nContractions = tContractions(n-a, arrays.array11, arrays.array22)
			if(len(pContractions) == 0):
				pContractions = [TypeContract()]
			if(len(nContractions) == 0):
				nContractions = [TypeContract()]
			for b in pContractions:
				for c in nContractions:
					contraction = Contract(b, c)
					contracts = contracts + [contraction]
		return contracts
	else:
		return []



def avoid_mc_1(conts, term1):         

    """ get rid of those contractions between active and inactive indices """
    """ 'conts' is  a list of Contract """

    result = []
    term = copy.deepcopy(term1)
    op1 = term.uOperators[0]
    op2 = term.uOperators[1]
    u1 = op1.upperIndicees
    u2 = op2.upperIndicees
    L2 = op2.lowerIndicees
    L1 = op1.lowerIndicees
    for cc in conts:
        h1 = cc.hTypeContract.indexList1
        h2 = cc.hTypeContract.indexList2
        p1 = cc.pTypeContract.indexList1
        p2 = cc.pTypeContract.indexList2
        flag = 0
        if(flag == 0):
            if(len(h1) != len(h2) or len(h1) > len(u1) or len(h1) > len(L2)):
                printterm(term1)
                print h1
                print h2
                pdb.set_trace()
                return
            if(len(p1) != len(p2) or len(p1) > len(u2) or len(p1) > len(L1)):
                printterm(term1)
                print p1
                print p2
                pdb.set_trace()
                return
        if(flag == 0):
            for xx in range(len(h1)):
                if(u1[h1[xx]].type == 'h' and u1[h1[xx]].active == -1 and L2[h2[xx]].type == 'm'):
                    flag = 1
                    break
                elif(L2[h2[xx]].type == 'h' and L2[h2[xx]].active == -1 and u1[h1[xx]].type == 'm'):
                    flag = 1
                    break
        if(flag == 0):
            for xx in range(len(p1)):
                if(L1[p1[xx]].type == 'p' and L1[p1[xx]].active == -1 and u2[p2[xx]].type == 'e'):
                    flag = 1
                    break
                elif(u2[p2[xx]].type == 'p' and u2[p2[xx]].active == -1 and L1[p1[xx]].type == 'e'):
                    flag = 1
                    break
        if(flag == 0):
                result += [cc]
    return result



def loopsToJoin(loops, a):

        """ tells the position of the loops that will be joined at a position 'a' """

	L = []
	for elem in loops:
		if(elem.startP == a or elem.endP == a):
			L = L + [loops.index(elem)]
		else:
			pass
	return L






def join(loop1, loop2):

        """ joins two loops """

	joint = Loop()
	joint.middles = loop1.middles + loop2.middles
	joint.type = loop2.type
        if(loop1.type is not loop2.type):
	        if(loop1.endP == loop2.endP):
			joint.middles = joint.middles + [loop1.endP]
			joint.startP = loop2.startP
			joint.endP = loop1.startP
		elif(loop1.startP == loop2.startP):
			joint.middles = joint.middles + [loop1.startP]
			joint.startP = loop1.endP
			joint.endP = loop2.endP

	elif(loop1.type == loop2.type):
		if(loop1.startP == loop2.endP):
			joint.middles = joint.middles + [loop1.startP]
			joint.startP = loop2.startP
			joint.endP = loop1.endP
		elif(loop1.endP == loop2.startP):
			joint.middles = joint.middles + [loop1.endP]
			joint.startP = loop1.startP
			joint.endP = loop2.endP

	return joint





def makeLoops(contraction, length):

        """ makes loops from 'contraction' which is an object of 'Contract' """
        """ 'length' is the number of pairs of labels in first operator """

	loops = []
	c11 = contraction.hTypeContract.indexList1
	c12 = contraction.pTypeContract.indexList2
	c21 = contraction.pTypeContract.indexList1
	c22 = contraction.hTypeContract.indexList2
	for a in range(len(c11)):
		loop = Loop("h", c11[a], c22[a]+length, [])
		loops = loops + [loop]
	for b in range(len(c21)):
		loop = Loop("p", c21[b], c12[b]+length, [])
		loops = loops + [loop]
	joints = []
	for c in c11:
		if(c in c21):
			joints = joints + [c]
	for d in c12:
		if(d in c22):
			joints = joints + [d+length]
	for e in joints:
		l = loopsToJoin(loops, e)
		if(len(l) < 2):
			pass
		else:
			joint = join(loops[l[0]], loops[l[1]])
			loops[l[0]] = joint
			loops = rmEl(loops[l[1]], loops)
	return loops






def coeff(contraction, len):

        """ gives the constant coefficient (number) for a contraction, which is determined by No. of loops and No. of hole labels """
        """ 'len' is the number of label pairs in first operator """

	l11 = contraction.hTypeContract.indexList1
	l12 = contraction.pTypeContract.indexList2
	l21 = contraction.pTypeContract.indexList1
	l22 = contraction.hTypeContract.indexList2
	coefficient = 1
        if(system_switch['use spin orbital'] == 0):
            Loop_coeff = -2
        elif(system_switch['use spin orbital'] == 1):
            Loop_coeff = -1
        else:
            print "\n\n   Unspecified value. coeff@real.py\n\n"
            pdb.set_trace()
	for a in l11:
		coefficient = coefficient*-1
	loops = makeLoops(contraction, len)
	for a in loops:
		if(a.startP == a.endP):
		    coefficient = coefficient*Loop_coeff
		else:
		    pass
	return coefficient






def generalIndex(ind1, ind2):          

        """ ind1 and ind2 are contracted; this function determine the right order to make the first index more general; """
        """ the key point is to pass all information in the 'general' index to the less general index """
        """ because in following steps the 'general' one will be replaced by the less general one """

        if ind1.type == ind2.type == "g":
            if ind1.fix == ind2.fix == 1: set_trace()
            elif ind1.fix == 1: return [ind2, ind1]
            elif ind2.fix == 1: return [ind1, ind2]
            else: return [ind1, ind2]
        elif ind1.type == "g":
            if(ind1.fix == 1): set_trace()
            else: return [ind1, ind2]
        elif ind2.type == "g":
            if ind2.fix == 1: set_trace()
            else: return [ind2, ind1]
        elif ind1.fix == 1:             # first take care of fix index -> 
            if(ind2.fix == 1):
                print "Problem in generalIndex()!!\n\n"; set_trace()
            elif dic_type_value[ind1.type] > dic_type_value[ind2.type]:
                return [ind2, ind1]
            elif dic_type_value[ind1.type] == dic_type_value[ind2.type] and ind1.active > ind2.active:
                ind2.fix = 1
                ind2.att = ind1.att
                return [ind1, ind2]
            elif dic_type_value[ind1.type] == dic_type_value[ind2.type] and ind1.active < ind2.active:
                return [ind2, ind1]
            elif dic_type_value[ind1.type] < dic_type_value[ind2.type]:
                ind2.att = ind1.att
                ind2.fix = 1
                return [ind1, ind2]  
            else:
                return [ind2, ind1]
        elif ind2.fix == 1:
            if dic_type_value[ind2.type] > dic_type_value[ind1.type]:
                return [ind1, ind2]
            elif dic_type_value[ind2.type] == dic_type_value[ind1.type] and ind2.active > ind1.active:
                ind1.att = ind2.att
                ind1.fix = 1
                return [ind2, ind1]
            elif dic_type_value[ind2.type] == dic_type_value[ind1.type] and ind2.active < ind1.active:
                return [ind1, ind2]
            elif dic_type_value[ind2.type] < dic_type_value[ind1.type]:
                ind1.att = ind2.att
                ind1.fix = 1
                return [ind2, ind1]
            else:
                return [ind1, ind2]
        elif ind1.type == "m":        # if there is no 'g' type of index, then 'm' must be special ->
            return [ind2, ind1]
        elif ind2.type == "m":
            return [ind1, ind2]    # if there is neither 'g' or 'm', 'e' must be speical ->
        elif ind1.type == "e":
            return [ind2, ind1]
        elif ind2.type == "e":
            return [ind1, ind2]    # if there is no fixed or 'g' or 'm' or 'e', then inactive must be special
        elif ind1.type == ind2.type == 'h':
            if ind1.active == -1 and ind2.active != -1: return [ind2,ind1]
            else: return [ind1, ind2]
        elif ind1.type == ind2.type == 'p':
             if ind1.active == -1 and ind2.active != -1: return [ind2, ind1]
             else: return [ind1, ind2]
        else:
            print "Problem in generalIndex()!!\n\n"
            set_trace()






def makeMat(contraction, term1):    

        """ replace 'general' contracted indices in matrixelement of a term by less general ones, outputs matelements """

	c11 = contraction.hTypeContract.indexList1
	c12 = contraction.pTypeContract.indexList2
	c21 = contraction.pTypeContract.indexList1
	c22 = contraction.hTypeContract.indexList2
        term = copy.deepcopy(term1)
	o1 = term.uOperators[0]
        o2 = term.uOperators[1]
	matElems = []
	gen = []
	holPar = []
	for a in range(len(c11)):
		i = generalIndex(o1.upperIndicees[c11[a]], o2.lowerIndicees[c22[a]])
		gen = gen + [i[0]]
		holPar = holPar + [i[1]]
	for a in range(len(c21)):
		i = generalIndex(o1.lowerIndicees[c21[a]], o2.upperIndicees[c12[a]])
		gen = gen + [i[0]]
		holPar = holPar + [i[1]]                                       # Now we have classified general and special indices which are from 'operators'. Next we will    
        term = reset(term, gen + holPar, gen + holPar)                         # compare with the indices in matelement to determine which index to keep. But it     
	for count in range(len(term.coefficient.matElement)):                  # turns out that the same index in operator and in mat is not pointing to 
	        uInd = term.coefficient.matElement[count].matUpperIndicees     # the same address, so we need to UPDATE the information of indices in mats, because by doing 
	        lInd = term.coefficient.matElement[count].matLowerIndicees     # generalIndex() some important information is coded. This is why we use reset() here.
		upperInd = []
		lowerInd = []
		for a in range(len(uInd)):
			if(uInd[a].inArray(gen)):
				upperInd = upperInd + [holPar[uInd[a].gvIndexNo(gen)]]
			else:
				upperInd = upperInd + [uInd[a]]
			if(lInd[a].inArray(gen)):
				lowerInd = lowerInd + [holPar[lInd[a].gvIndexNo(gen)]]
			else:
				lowerInd = lowerInd + [lInd[a]]
		matElems = matElems + [MatElement(term.coefficient.matElement[count].name, upperInd, lowerInd)]
	return matElems



def makeTerm(contraction, term):      

        """ makes the whole term for a particular contraction: replace 'general' by less general contracted indices """
    
	c11 = contraction.hTypeContract.indexList1
	c12 = contraction.pTypeContract.indexList2
	c21 = contraction.pTypeContract.indexList1
	c22 = contraction.hTypeContract.indexList2
	uppbk0 = []
	lowbk0 = []
	uppbk1 = []
	lowbk1 = []
        o1 = copy.deepcopy(term.uOperators[0])
        o2 = copy.deepcopy(term.uOperators[1])
	mats = makeMat(contraction, term)
	coef = term.coefficient.const
	coef = coef*coeff(contraction, len(o1.upperIndicees))
	loops = makeLoops(contraction, len(o1.upperIndicees))
	uInd = o1.upperIndicees + o2.upperIndicees
	lInd = o1.lowerIndicees + o2.lowerIndicees
	uppers = []
	lowers = []
	middle = []
	for a in loops:
		middle = middle + a.middles
		if(a.startP == a.endP):
			middle = middle + [a.endP]
		elif(a.type == "p"):
			uppers = uppers + [a.endP]
			lowers = lowers + [a.startP]	
		elif(a.type == "h"):
			uppers = uppers + [a.startP]
			lowers = lowers + [a.endP]
	upperInd = []
	lowerInd = []
        delete = []
	for b in range(len(uInd)):
		if(b in middle): pass
		else:
			if(b in delete): pass
			elif(b in uppers):
			    if(b in lowers): pass
			    else:
			        lowerInd = lowerInd + [lInd[b]]
			        upperInd = upperInd + [uInd[lowers[uppers.index(b)]]]
			        delete = delete + [lowers[uppers.index(b)]]
			elif b in lowers:
				upperInd = upperInd + [uInd[b]]
				lowerInd = lowerInd + [lInd[uppers[lowers.index(b)]]]
				delete = delete + [uppers[lowers.index(b)]]
                        else:
		            upperInd = upperInd + [uInd[b]]
		            lowerInd = lowerInd + [lInd[b]]
	return Term(Coefficient(coef, mats), [Uoperator(upperInd, lowerInd)])





#**************************************************
#                      show                       *
#**************************************************




def makeIndTerm(term1):

        """ makes an Indicees object for the whole term """

	ind = []
        term = copy.deepcopy(term1)
	for a in term.uOperators:
		ind = ind + [makeIndexSet(a.upperIndicees, a.lowerIndicees)]
	for b in term.coefficient.matElement:
		ind = ind + [makeIndexSet(b.matUpperIndicees, b.matLowerIndicees)]
	return Indicees(ind)





def arrangeIndicees(ind, term):

        """ arranges the labels from 'ind' object to an object of Term (class), first operators then matrixelements """

	for n in range(len(ind.indexSets)):
		set = ind.indexSets[n]
                if(n < len(term.uOperators)):
			for a in set.indexPairs:
				term.uOperators[n].upperIndicees[set.indexPairs.index(a)] = a.upperIndex
				term.uOperators[n].lowerIndicees[set.indexPairs.index(a)] = a.lowerIndex
		else:
			for a in set.indexPairs:
				term.coefficient.matElement[n-len(term.uOperators)].matUpperIndicees[set.indexPairs.index(a)] = a.upperIndex
				term.coefficient.matElement[n-len(term.uOperators)].matLowerIndicees[set.indexPairs.index(a)] = a.lowerIndex
	return term



def solveforncontrac(term1, n):

        """ solves the first two operators of the term for 'n' contractions """
        """ term1 has two operators in term1.uOperators """

        result = []
        term = copy.deepcopy(term1)
	if(n == 0):
		opinter = Uoperator(term.uOperators[0].upperIndicees + term.uOperators[1].upperIndicees,term.uOperators[0].lowerIndicees + term.uOperators[1].lowerIndicees)
                co = Coefficient(term.coefficient.const, term.coefficient.matElement)
                inter = Term(co, [opinter])
              	return [inter]
	else:
		conts = contractions(n, makeContableArrays(term))
                conts = avoid_mc_1(conts, term)
		if(len(conts) == 0):
                        #print "no contraction formed, return emply list \n"
			return []
		else:
			i = 1
			for z in conts:
				temp = makeTerm(z, term)
                                if(n == 1 and 0):
                                    printterm(term1)
                                    printterm(temp)
                                    pdb.set_trace()
				result = result + [temp]
				i = i + 1
			return result



def maxContraction(term):              

    """ tells the maximum number of contractions possible in a product of two operators """
  
    arr = makeContableArrays(term)
    return min(len(arr.array11), len(arr.array22)) + min(len(arr.array12), len(arr.array21))



def solvetwoop(term):
 
    """ solves the term for two operators; essence: expand A.B in normal order according to GWT, assume A and B themselves already in normal order form """
 
    if(system_switch['normal order'] == 'MN'):
        return  mn_gwt_solvetwoop(term)
    elif(system_switch['normal order'] == 'MN2'):
        return  mn2_gwt_solvetwoop(term)
    elif(system_switch['normal order'] == 'KM'):
        return  km_gwt_solvetwoop(term)
    elif(system_switch['normal order'] == 'T'):
        return  t_gwt_solvetwoop(term)
    elif(system_switch['normal order'] == 'V'):
        return  v_gwt_solvetwoop(term)
    else:
        print "\n  normal order undefined"
        pdb.set_trace()



def t_gwt_solvetwoop(term):
    result = []
    if(len(term.uOperators) == 1 ):
        print "\n this term has only one operator, no contraction formed \n"
        return [term]
    else:
        for a in range(maxContraction(term)+1):
            result = result + solveforncontrac(term, a)
        return result




def v_gwt_solvetwoop(term):
    """ this is implemented from tweaking km_gwt_solvetwoop """
    if(system_switch['explicit spin'] != 0):
        print "\n Genuine vacuum normal order: 'explicit spin' should be set to 0\n"
        print system_switch['explicit spin']
        set_trace()
    elif(system_switch['GWT CU truncate rank'] != 1):
        print "\n Genuine vacuum normal order: 'GWT CU truncate rank' should be set to 1\n"
        print system_switch['GWT CU truncate rank']
        set_trace()
    elif(system_switch['Allow Internal Contraction'] != 0):
        print "\n Genuine vacuum normal order: 'Allow Internal Contraction' should be set to 0\n"
        print system_switch['Allow Internal Contraction']
        set_trace()
    else:
        res = km_gwt_solvetwoop(term)
        res = km_gwt_decode_eta_s(res)
        res = km_gwt_decode_delta_s(res)
        if(system_switch['use spin orbital'] == 0):
            res = select_terms_rank_n(res, [], ['Lambda1'])
        else:
            res = select_terms_rank_n(res, [], ['lambda1'])
        return res














#--------------------------------------------------------------------------------------------------------------------------------------
                
#                                               KM-GWT
                
#--------------------------------------------------------------------------------------------------------------------------------------
                






class KM_GWT_Separator:
    def __init__(self, start = "Y"):
        self.start = type








# This section is for Kutzelnigg & Mukherjee normal order theory

                
def km_gwt_solvetwoop(term):
    """                
        solves the term for two operators; essence: expand A.B in normal order according to GWT, assume A and B themselves already in normal order form 
    """
    debug = 0
    if(debug):
        print " \n\n debugging km_gwt_solvetwoop:\n\n"
    result = []
    if(len(term.uOperators) == 1 ):
        print "\n this term has only one operator, no contraction formed \n"
        return [term]
    else:
        maxcon = len(term.uOperators[0].upperIndicees + term.uOperators[1].upperIndicees)    # this is the maximal number of contraction allowed 
        if(debug):
            print " look at 2 contraction to debug:"
            a = 2
            pdb.set_trace()
            result +=  km_gwt_solveforncontrac(term, a)
        else:
            for a in range(maxcon + 1):
#                if(a == 3):
#                    pdb.set_trace()
                result += km_gwt_solveforncontrac(term, a)
        return result




def km_gwt_solveforncontrac(term1, n):
    """
        solves the first two operators of the term for 'n' contractions 
        term1 has two operators in term1.uOperators 
    """
    debug = 0
    if(debug):
        print "\n debug km_gwt_solveforncontrac \n"
    result = []
    term = copy.deepcopy(term1)
    LL = len(term.uOperators[0].upperIndicees)
    opinter = Uoperator(term.uOperators[0].upperIndicees + term.uOperators[1].upperIndicees,term.uOperators[0].lowerIndicees + term.uOperators[1].lowerIndicees)
    co = Coefficient(term.coefficient.const, term.coefficient.matElement)
    if(n == 0):
        inter = Term(co, [opinter])
        return [inter]
    else:
        if(debug):
            print "\n print out n and LL\n"
            print n
            print LL
            print "\n"
            pdb.set_trace()
        conts = km_gwt_contractions(n, LL, opinter)    # for km_gwt, contractin doesn't depend on the matrixelements and coeff,
        if(debug and n == 2):                                    
            #pdb.set_trace()
            for XYZ in conts:
                print_num_operator_s(XYZ) 
                #pdb.set_trace()
            pdb.set_trace()
        if(len(conts) == 0):                           # so we treat contractions SEPARATELY and add the orginial coeff at the end after contractions done
            #print "no contraction formed, return emply list \n"
            return []
        else:
            i = 1
            if(0 and debug):
                print "\n printout conts:"  
                for z in conts:
                    print_num_operator_s(z)
                print "\n"
            if(debug):
                print "\n number of elements in conts: " + str(len(conts)) + "\n"
            for z in conts:
                if(debug):
                    print "\n printout z:"
                    print_num_operator_s(z)
                    print "\n"
                if(system_switch['use spin orbital'] == 0):
                    temp = km_gwt_sf_makeTerm(z, opinter)
                elif(system_switch['use spin orbital'] == 1):
                    temp = km_gwt_so_makeTerm(z, opinter)
                else:
                    print "\n\n in km_gwt_solveforncontrac; system_switch['use spin orbital'] undefined \n\n"
                    pdb.set_trace()
                if(debug):
                    print "\n print out temp:"
                    printterms(temp)
                    pdb.set_trace()
                result = result + temp
            for y in result:
                y.coefficient.const *= co.const
                y.coefficient.matElement += co.matElement
            return result



def km_gwt_contractions(n, LL, opinter):
    """
     find all possible n contractions, return a list of objects of the class Uoperator, whose upper and lower indicees are the positions of indices contracted
     LL is the number of (upper) indices from the first operator 
     for both SR theory and MR theory, if we know which indices are contracted, the other thing to determine is how to partion them
     for one particular partition pattern, there are different possible combinations
     for SR theory, the partion is unique: [1, 1, 1, ...], for MR theory, multiple contraciton is possible, so partition pattern are  versatile
     the special thing about spin free theory: we maximize the matching between upper and lower indicees
     but we can indeed use the spin-free contraction rule for spin orbital theory at this step, just a bit more work is done;
     on the other hand, we can fix the partition pattern to be  [1, 1, 1, ...] and then check the algorithm
    """
    debug = 0
    if(debug):
        print "\n\n debug km_gwt_contractions:\n"
    check_with_sr = 0
    result = []
    valid_upp = valid_low = range(len(opinter.upperIndicees))
    upp_contrac_list = combi(n, valid_upp)                 # all contracted indices
    low_contrac_list = combi(n, valid_low)
    Parts = group_partition_int(n)
    parts = []
    if(check_with_sr):                                     # check with SR theory
        A = []
        for xxa in range(n):
            A += [1]                                       # A = [1, 1, 1, ...]
        parts = [ Uoperator([1],[n]) ]
    elif(debug and 0):
        print "\n n = 3, set parts = [ [1, 2] ] to debug\n"
        parts = [ Uoperator([1,2],[1,1]) ]  
    elif(system_switch['GWT CU truncate rank'] > 0):                     # truncate cumulants
        for xy in Parts:
            ok = 1
            for yx in xy.upperIndicees:
                if(yx > system_switch['GWT CU truncate rank']):
                    ok = 0
                    break
            if(ok == 1):
                parts += [xy]
    else:
        parts = Parts
    if(debug):
        print "\n parts:\n"
        print_part_s(parts)
        print "\n"
        pdb.set_trace()
    for upp_con in upp_contrac_list:
        for low_con in low_contrac_list:                   # remove those contractions that all contracted indices are from one operator
            if(not both_in_one_op(upp_con, low_con, LL) or system_switch['Allow Internal Contraction']):   
                for part in parts:                         # deal with one particular partition pattern
                    result_component = km_gwt_contractions_partition(LL, upp_con, low_con, part)
                    if(0  and n == 3):
                        print "part:"
                        print_num_operator(part)
                        print_num_operator_s_s(result_component)
                        pdb.set_trace()
                    result += result_component
                    #if(debug and upp_con == [1, 2, 3] and low_con == [0, 1,2] and part.upperIndicees == [1,2] and part.lowerIndicees == [1,1] ):
                    if(debug):
                        print "\n part, upp_con, low_con:\n"
                        print_part(part)
                        print upp_con
                        print low_con
                        print "\n" 
                        for dd in result_component:
                            print_num_operator_s(dd)
                        print " \n done \n"
                        pdb.set_trace()
    return result


def print_part_s(parts):
    for xx in parts:
        print_part(xx)
    return



def print_part(part):
    strr = ""
    for xx in part.upperIndicees:
        for yy in range(part.lowerIndicees[0]):
            strr += str(xx) + ' '
    print strr



def km_gwt_contractions_partition(LL, upp_con, low_con, part):
    """ 
        give all possible contraction patterns for a particular partion pattern and known contracted upper and lower indicees 
        return e.g., (if upp_con = [0,1,3] & low_con = [2,4,5], part = op(1,2; 1,1)),  [ [ Uop(3;4), Uop(0,1; 2,5)], ..., ...  ] 
        maximize matching 
    """
    result = []
    debug = 0
    if(debug):
        print "\n debug km_gwt_contractions_partition: \n"
        pdb.set_trace()
    if(len(part.upperIndicees) == 0):
        result = []
    elif(len(part.upperIndicees) == 1):
        result = km_gwt_contractions_partition_section(LL, upp_con, low_con, part.upperIndicees[0], part.lowerIndicees[0])         # RECURSION start
    else:
        result_part_1 = km_gwt_contractions_partition_section(LL, upp_con, low_con, part.upperIndicees[0], part.lowerIndicees[0])  
        if(debug):
            for xyz in result_part_1:
                print_num_operator_s(xyz)
            pdb.set_trace()
        for ops in result_part_1:                                              # e.g. ops = [Uop(3;4), Uop(1;2), Uop(5, 3)]
            upp_section = []
            low_section = []
            for xx in ops:
                upp_section += xx.upperIndicees
                low_section += xx.lowerIndicees
            upp_delete_offsets = []
            low_delete_offsets = []
            for xxx in upp_section:
                upp_delete_offsets += [upp_con.index(xxx)]
            for yyy in low_section:
                low_delete_offsets += [low_con.index(yyy)] 
            upp_con_left = deletes(upp_con, upp_delete_offsets)
            low_con_left = deletes(low_con, low_delete_offsets)                # RECURSION point
            contrac_from_other_sec = km_gwt_contractions_partition( LL, upp_con_left, low_con_left, Uoperator(part.upperIndicees[1:], part.lowerIndicees[1:]) )
            for zzz in contrac_from_other_sec:
                result += [ ops + zzz ] 
    if(debug):
        pdb.set_trace()
        for xyz in result:
            print_num_operator_s(xyz)
    return result


def km_gwt_contractions_partition_section(LL, upp_con, low_con, section_size, section_num):
    if(system_switch['use spin orbital'] == 0):
        return km_gwt_contractions_partition_section_sf(LL, upp_con, low_con, section_size, section_num)
    elif(system_switch['use spin orbital'] == 1):
        return km_gwt_contractions_partition_section_so(LL, upp_con, low_con, section_size, section_num)
    else:
        print "\n global parameter undefined"
        pdb.set_trace()



def km_gwt_contractions_partition_section_so(LL, upp_con, low_con, section_size, section_num):
    debug = 0
    if(debug):
        print "\n debug km_gwt_contractions_partition_section \n"
    result = []
    upp_possible_combi = unique_section_divide(upp_con, section_size, section_num) 
    low_possible_combi = permu_unique_section_divide(low_con, section_size, section_num)
    if(debug):
        pdb.set_trace()
        print upp_possible_combi
        print low_possible_combi   
    for xx in upp_possible_combi:
        for yy in low_possible_combi:
            ops = []
            self_contraction = 0
            for zz in range(section_num):
                if(system_switch['Allow Internal Contraction'] == 0):
                    if(both_in_one_op( xx[zz], yy[zz], LL)):                                           
                        self_contraction = 1                                                           # get rid of self_contraction
                        break
                    else:
                        ops += [ Uoperator( xx[zz], max_match(xx[zz], yy[zz]) ) ]                      # maximize matching
                else:  
                    ops += [ Uoperator( xx[zz], max_match(xx[zz], yy[zz]) ) ]                          # from every pair(xx,yy), we get a 'ops'
            if(self_contraction == 1 and system_switch['Allow Internal Contraction'] == 0):
                pass
            else:
                result +=  [ copy.deepcopy(ops) ]
    return result



def km_gwt_contractions_partition_section_sf(LL, upp_con, low_con, section_size, section_num):
    debug = 0
    if(debug):
        print "\n debug km_gwt_contractions_partition_section \n"
    result = []
    upp_possible_combi = unique_section_divide(upp_con, section_size, section_num)
    if(section_size == 1):
        low_possible_combi = permu_unique_section_divide(low_con, section_size, section_num)
    else:
        low_possible_combi = inter_out_permu_unique_section_divide(low_con, section_size, section_num)# outside permu followed by internal permu 
    if(debug):
        print section_size
        print section_num
        print upp_possible_combi
        print low_possible_combi
        pdb.set_trace()
    for xx in upp_possible_combi:
        for yy in low_possible_combi:
            ops = []
            ok = 1
            for zz in range(section_num):
                xyx = max_match(xx[zz], yy[zz])
                if(xyx != yy[zz]):                                                                     # make sure common indices (upper&lower) are lined up
                    ok = 0
                    break
                elif(system_switch['Allow Internal Contraction'] == 0):
                    if(both_in_one_op( xx[zz], yy[zz], LL)):
                        ok = 0                                                                         # get rid of self_contraction
                        break
                    else:
                        ops += [ Uoperator( xx[zz], xyx ) ]                                            # maximize matching
                else:
                    ops += [ Uoperator( xx[zz], xyx ) ]                                                # from every pair(xx,yy), we get a 'ops'
            if(ok == 1):
                result +=  [ copy.deepcopy(ops) ]
    return result





def km_gwt_sf_makeTerm(con_pattern, opinter):
    """
     from contraction pattern make the new terms: contractions become eta/Gamma/Lambda 
     con_pattern contains the positions of indices contracted, LL is the number of (upper) indices from the first op, opinter is the merged op 
     task: determine the number of km_gwt_contractions_partition_section(LL, oops -> each loop coresponds to a (-2) multiplication factor)
           determine the correspondence of uncontracted indices -> determine the final operator
           chagnge the contraction to corresponding matrix elements
    """
    debug = 0
    if(debug):
        print "\n\n debugging km_gwt_sf_makeTerm \n\n"
        pdb.set_trace()
    valid = 1
    loop_coeff = -2
    result = Term(Coefficient(1.0, []), [])
    total_L = len(opinter.upperIndicees)
    upp_con = []
    low_con = []
    for xx in con_pattern:
        upp_con += xx.upperIndicees
        low_con += xx.lowerIndicees
    uncontrac_inds_pos_pair = pair_uncontra(upp_con, low_con, total_L)
    final_op = Uoperator([], [])
    for xx in range(len(uncontrac_inds_pos_pair[0])):
        final_op.upperIndicees += [ opinter.upperIndicees[ uncontrac_inds_pos_pair[0][xx] ] ]
        final_op.lowerIndicees += [ opinter.lowerIndicees[ uncontrac_inds_pos_pair[1][xx] ] ]
    result.uOperators = [final_op]                                             # get operator
    num_loops = find_num_loops(upp_con, low_con)
    if(debug):
        pdb.set_trace()
        print " \n uncontrac_inds_pos_pair:"
        print uncontrac_inds_pos_pair[0]
        print uncontrac_inds_pos_pair[1]
        print "\n print opinter:"
        printOperator(opinter)
        print "\n final_op"
        printOperator(final_op)
        print num_loops
        pdb.set_trace()
    for A in range(num_loops):
        result.coefficient.const *= loop_coeff                                 # coeff contribution from loops
    for yy in con_pattern:                                                     # contractions contribution
        mat_yy = km_gwt_sf_make_mat_one(yy, opinter)
        if(mat_yy[0] == 1):
            result.coefficient.const *= mat_yy[1]
            result.coefficient.matElement += [mat_yy[2]]                         
        elif(mat_yy[0] == 0):
            valid = 0
            break
        else:
            pdb.set_trace()
    if(valid == 1):
        return [result]
    elif(valid == 0):
        return []
    else:
        pdb.set_trace()
        return 


def km_gwt_so_makeTerm(con_pattern, opinter):
    """
      from contraction pattern make the new terms: contractions become eta/Gamma/Lambda 
      con_pattern contains the positions of indices contracted, LL is the number of (upper) indices from the first op, opinter is the merged op 
      task: determine the number of loops -> each loop coresponds to a (-2) multiplication factor
            determine the correspondence of uncontracted indices -> determine the final operator
            chagnge the contraction to corresponding matrix elements
    """
    debug = 0
    if(debug):
        print "\n\n debugging km_gwt_so_makeTerm \n\n"
        pdb.set_trace()
    valid = 1
    loop_coeff = -1
    result = Term(Coefficient(1.0, []), [])
    total_L = len(opinter.upperIndicees)
    upp_con = []
    low_con = []
    for xx in con_pattern:
        upp_con += xx.upperIndicees
        low_con += xx.lowerIndicees
    uncontrac_inds_pos_pair = pair_uncontra(upp_con, low_con, total_L)
    final_op = Uoperator([], [])
    for xx in range(len(uncontrac_inds_pos_pair[0])):
        final_op.upperIndicees += [ opinter.upperIndicees[uncontrac_inds_pos_pair[0][xx]] ]
        final_op.lowerIndicees += [ opinter.lowerIndicees[uncontrac_inds_pos_pair[1][xx]] ]
    result.uOperators = [final_op]                                             # get operator
    num_loops = find_num_loops(upp_con, low_con)
    for A in range(num_loops):
        result.coefficient.const *= loop_coeff                                 # coeff contribution from loops
    for yy in con_pattern:                                                     # contractions contribution
        mat_yy = km_gwt_so_make_mat_one(yy, opinter)
        if(mat_yy[0] == 1):
            result.coefficient.const *= mat_yy[1]
            result.coefficient.matElement += [mat_yy[2]]                           
        elif(mat_yy[0] == 0):
            valid = 0
            break
        else:
            pdb.set_trace()
    if(valid == 1):
        return [result]
    elif(valid == 0):
        return []
    else:
        pdb.set_trace()
        return


def km_gwt_sf_make_mat_one(one_con, opinter):
    if(system_switch['Allow Internal Contraction'] == 0):
        return km_gwt_sf_make_mat_one_no_inter_contraction(one_con, opinter)
    elif(system_switch['Allow Internal Contraction'] == 1):
        return km_gwt_sf_make_mat_one_with_inter_contraction(one_con, opinter)
    else:
        print "\n condition undefined \n"
        pdb.set_trace()
        return


def km_gwt_sf_make_mat_one_no_inter_contraction(one_con, opinter):
    """
        the argument one_con is a particular one/multiple contraction, e.g. one_con = op(1,2; 3,4), return [const, mat] 
    """ 
    debug = 0
    valid = 1                                              # when system_switch['EXCLUDE Particle Ind from CU'] == 1 and there is 'p' label in cumulants or DM's, valid = 0        
    const = 1.0
    eta_coeff = 0.5
    mis_match_coeff = -0.5                                 # '-1' takes care of sign, '0.5' is because here we use SPATIAL DM, while one contraction corresponds
    upp_inds = []                                          # to a spin DM 
    low_inds = []
    for xx in range(len(one_con.upperIndicees)):
        upp_inds += [opinter.upperIndicees[ one_con.upperIndicees[xx] ]]
        low_inds += [opinter.lowerIndicees[ one_con.lowerIndicees[xx] ]]
    if(len(one_con.upperIndicees) == 0):
        print "\n\n in km_gwt_make_mat_one, point 1\n\n"
        pdb.set_trace()
    elif(len(one_con.upperIndicees) == 1):
        if(one_con.upperIndicees[0] < one_con.lowerIndicees[0]):
            const *= mis_match_coeff                          
            mat = MatElement('Lambda1', upp_inds, low_inds)
            if(system_switch['EXCLUDE Particle Ind from CU'] == 1):              # get rid of cumulants containing particle labels
                types = counttype(upp_inds + low_inds)
                if(types[0] + types[1] > 0):
                    valid = 0
        elif(one_con.upperIndicees[0] > one_con.lowerIndicees[0]):
            const *= eta_coeff
            mat = MatElement('eta', upp_inds, low_inds)
            if(system_switch['EXCLUDE Particle Ind from CU'] == 1):              # in this case, \eta^a_i == 0 
                if(upp_inds[0].type in ['p', 'P', 'e', 'E', 'a'] and low_inds[0].type in ['h', 'H', 'm', 'M', 'i']):
                    valid = 0
                elif(upp_inds[0].type in ['h', 'H', 'm', 'M', 'i'] and low_inds[0].type in ['p', 'P', 'e', 'E', 'a']): 
                    valid = 0
        else:
            print "\n\n in km_gwt_sf_make_mat_one, point 2 \n\n"
            pdb.set_trace()
    else:
        mat = MatElement('Lambda', upp_inds, low_inds)
        for xx in range(len(one_con.upperIndicees)):
            if(one_con.upperIndicees[xx] != one_con.lowerIndicees[xx]):
                const *= mis_match_coeff
        if(system_switch['EXCLUDE Particle Ind from CU'] == 1):                  # get rid of cumulants containing particle labels
            types = counttype(upp_inds + low_inds)
            if(types[0] + types[1] > 0):
                valid = 0
    return [valid, const, mat]



def km_gwt_sf_make_mat_one_with_inter_contraction(one_con, opinter):
    """
        I think this function is right for the spin-orbital case, may not be so for spin-free case
        one_con is a particular one/multiple contraction, e.g. one_con = op(1,2; 3,4), return [const, mat] 
        This function is written solely for the purpose of generating the relations between old and new ops/DM's/cumulants 
        so (1) internal contraction is allowd (2) for single contractions, all give a 'Lambda' contribution; mismatch contributes 
        -0.5, match contributes 1.0 to the final term; (3) for multiple contractions, the upper and lower indices in 'Lambda' have 
        been aligned up to the maximal extent, then every mismatch contributes -0.5 
    """
    print "\n this function could be wrong; check before using it \n"
    pdb.set_trace()
    debug = 0
    valid = 1                                              # when system_switch['EXCLUDE Particle Ind from CU'] == 1 and there is 'p' label in cumulants or DM's, valid = 0
    const = 1.0                                            # when internal contractions are allowed, we are expanding E^.._.., instead of the product of two E operators 
    mis_match_coeff = -0.5                                 # so there should be no 'eta'
    match_coeff = 1.0
    upp_inds = []
    low_inds = []
    for xx in range(len(one_con.upperIndicees)):
        upp_inds += [opinter.upperIndicees[ one_con.upperIndicees[xx] ]]
        low_inds += [opinter.lowerIndicees[ one_con.lowerIndicees[xx] ]]
    if(len(one_con.upperIndicees) == 0):
        print "\n\n in km_gwt_make_mat_one, point 1\n\n"
        pdb.set_trace()
    elif(len(one_con.upperIndicees) == 1):
        mat = MatElement('Lambda1', upp_inds, low_inds)
        if(system_switch['EXCLUDE Particle Ind from CU'] == 1):                  # get rid of cumulants containing particle labels
            types = counttype(upp_inds + low_inds)
            if(types[0] + types[1] > 0):
                valid = 0
        if(one_con.upperIndicees[0] != one_con.lowerIndicees[0]):
            const = mis_match_coeff
    else:
        mat = MatElement('Lambda', upp_inds, low_inds)
        for xx in range(len(one_con.upperIndicees)):
            if(one_con.upperIndicees[xx] != one_con.lowerIndicees[xx]):
                const *= mis_match_coeff
        if(system_switch['EXCLUDE Particle Ind from CU'] == 1):                  # get rid of cumulants containing particle labels
            types = counttype(upp_inds + low_inds)
            if(types[0] + types[1] > 0):
                valid = 0
    return [valid, const, mat]




def km_gwt_so_make_mat_one(one_con, opinter):
    """
        one_con is a particular one/multiple contraction, e.g. one_con = op(1,2; 3,4), return [const, mat] 
    """
    if system_switch['Allow Internal Contraction'] == 0:
        return km_gwt_so_make_mat_one_no_inter_contraction(one_con, opinter)
    elif system_switch['Allow Internal Contraction'] == 1:
        return km_gwt_so_make_mat_one_with_inter_contraction(one_con, opinter)
    else:
        print "\n condition undefined in km_gwt_so_make_mat_one \n"
        set_trace()
        raise Exception('condition undefined in km_gwt_so_make_mat_one@eq_lib')


def km_gwt_so_make_mat_one_no_inter_contraction(one_con, opinter):
    debug = 0
    valid = 1                                              # when system_switch['EXCLUDE Particle Ind from CU'] == 1 and there is 'p' label in cumulants or DM's, valid = 0
    const = 1.0
    eta_coeff = 1.0
    mis_match_coeff = -1
    upp_inds = []
    low_inds = []
    for xx in range(len(one_con.upperIndicees)):
        upp_inds += [opinter.upperIndicees[ one_con.upperIndicees[xx] ]]
        low_inds += [opinter.lowerIndicees[ one_con.lowerIndicees[xx] ]]
    if(len(one_con.upperIndicees) == 0):
        print "\n\n in km_gwt_make_mat_one \n\n"
        pdb.set_trace()
    elif(len(one_con.upperIndicees) == 1):
        if(one_con.upperIndicees[0] < one_con.lowerIndicees[0]):
            const *= mis_match_coeff
            mat = MatElement('lambda1', upp_inds, low_inds)
            if(system_switch['EXCLUDE Particle Ind from CU'] == 1):              # get rid of cumulants containing particle labels
                types = counttype(upp_inds + low_inds)
                if(types[0] + types[1] > 0):
                    valid = 0
        elif(one_con.upperIndicees[0] > one_con.lowerIndicees[0]):
            const *= eta_coeff
            mat = MatElement('eta', upp_inds, low_inds)
            if(system_switch['EXCLUDE Particle Ind from CU'] == 1):              # in this case, \eta^a_i == 0
                if(upp_inds[0].type in ['p', 'P', 'e', 'E', 'a'] and low_inds[0].type in ['h', 'H', 'm', 'M', 'i']):
                    valid = 0
                elif(upp_inds[0].type in ['h', 'H', 'm', 'M', 'i'] and low_inds[0].type in ['p', 'P', 'e', 'E', 'a']):
                    valid = 0
        else:
            print "\n\n in km_gwt_so_make_mat_one \n\n"
            set_trace()
            raise Exception("\n\n in km_gwt_so_make_mat_one \n\n")
    else:
        mat = MatElement('lambda', upp_inds, low_inds)
        for xx in range(len(one_con.upperIndicees)):
            if(one_con.upperIndicees[xx] != one_con.lowerIndicees[xx]):
                const *= mis_match_coeff
        if(system_switch['EXCLUDE Particle Ind from CU'] == 1):                  # get rid of cumulants containing particle labels
            types = counttype(upp_inds + low_inds)
            if(types[0] + types[1] > 0):
                valid = 0
    return [valid, const, mat]




def km_gwt_so_make_mat_one_with_inter_contraction(one_con, opinter):
    print "\n this function is not tested yet; but probably correct  \n"
    pdb.set_trace()
    debug = 0
    valid = 1                                              # when system_switch['EXCLUDE Particle Ind from CU'] == 1 and there is 'p' label in cumulants or DM's, valid = 0
    const = 1.0
    eta_coeff = 1.0
    mis_match_coeff = -1
    upp_inds = []
    low_inds = []
    for xx in range(len(one_con.upperIndicees)):
        upp_inds += [opinter.upperIndicees[ one_con.upperIndicees[xx] ]]
        low_inds += [opinter.lowerIndicees[ one_con.lowerIndicees[xx] ]]
    if(len(one_con.upperIndicees) == 0):
        print "\n\n in km_gwt_make_mat_one \n\n"
        pdb.set_trace()
    elif(len(one_con.upperIndicees) == 1):
        mat = MatElement('lambda1', upp_inds, low_inds)
        if(system_switch['EXCLUDE Particle Ind from CU'] == 1):                  # get rid of cumulants containing particle labels
                types = counttype(upp_inds + low_inds)
                if(types[0] + types[1] > 0):
                    valid = 0
        if(one_con.upperIndicees[0] != one_con.lowerIndicees[0]):
            const *= mis_match_coeff
    else:
        mat = MatElement('lambda', upp_inds, low_inds)
        for xx in range(len(one_con.upperIndicees)):
            if(one_con.upperIndicees[xx] != one_con.lowerIndicees[xx]):
                const *= mis_match_coeff
        if(system_switch['EXCLUDE Particle Ind from CU'] == 1):                  # get rid of cumulants containing particle labels
            types = counttype(upp_inds + low_inds)
            if(types[0] + types[1] > 0):
                valid = 0
    return [valid, const, mat]





def pair_uncontra(upp_con, low_con, total_L):

    """ return uncontracted ind pos [ [1,3], [2,4] ] """

    upp_left = []
    low_left = []
    for xx in range(total_L):
        if(xx in upp_con):
            pass
        else:
            upp_left += [xx]
            if(xx not in low_con):
                low_left += [xx]
            else:
                yy = upp_con[ low_con.index(xx) ]
                while(yy in low_con):
                    yy = upp_con[ low_con.index(yy) ]
                low_left += [yy]
    return [upp_left, low_left]






def find_num_loops(upp_con, low_con):
    result = 0
    counted = []
    for xx in range(len(upp_con)):
        upp = upp_con[xx]
        low = low_con[xx]
        if(upp in counted):
            pass
        elif(upp == low):              # this case takes place in multiple contractions or internal-conractions when expanding a V-no operator into a KM-no operator; not counted as a loop
            counted += [upp]
            pass
        else:
            while(low != upp and low in upp_con):
                if(low not in counted):
                    counted += [low]
                low = low_con[upp_con.index(low)]
            if(low == upp):
                result += 1
    return result




#****************************************************
#  replace 'eta' by 'delta' and 'gamma' or 'Gamma'  *
#****************************************************




def km_gwt_decode_contraction_s(terms):                                        # expand \eta in terms of delta and gamma/Gamma
    result = []                                                                # it is done automatically to get rid of terms including particle-index-containing cumulants
    for xx in terms:
        result += km_gwt_decode_contraction(xx)
    return result



def remove_p_cumulants(terms):                                                 # remove cumulants containing particle indices
    result = []
    for xx in terms:
        flag = 0
        for yy in xx.coefficient.matElement:
            if(yy.name in ['lambda', 'lambda1', 'Lambda', 'Lambda1', 'gamma', 'Gamma']):
                for zz in yy.matUpperIndicees + yy.matLowerIndicees:
                    if(zz.type in ['p', 'P', 'a', 'e', 'E']):
                        flag = 1
                        break
                if(flag == 1):
                    break
        if(flag == 0):
            result += [xx]
    return result



def remove_core_cumulants(terms):                                              # remove cumulants containing core indices
    result = []                                                                # convert 1-particle core DM to 1 or 2 
    for xx in terms:
        flag = 1
        for yy in xx.coefficient.matElement:
            if(yy.name in ['lambda', 'lambda1', 'Lambda', 'Lambda1', 'gamma', 'Gamma']):
                if(len(yy.matUpperIndicees) > 1 and yy.name in ['lambda', 'lambda1', 'Lambda', 'Lambda1']):
                    for zz in yy.matUpperIndicees + yy.matLowerIndicees:
                        if(zz.active == -1):
                            flag = 0
                            break
                elif(len(yy.matUpperIndicees) == 1):
                    if(yy.matUpperIndicees[0].active == -1 and yy.matLowerIndicees[0].active == 1):
                        flag = 0
                        break
                    elif(yy.matUpperIndicees[0].active == 1 and yy.matLowerIndicees[0].active == -1):
                        flag = 0
                        break        
                    elif(yy.matUpperIndicees[0].active == -1 or yy.matLowerIndicees[0].active == -1):
                        if(system_switch['use spin orbital'] == 1):
                            yy.name = 'delta'
                        elif(system_switch['use spin orbital'] == 1):
                            yy.name = 'delta'
                            xx.coefficient.const *= 2.0  
            if(flag == 0):
                break
        if(flag == 1):
            result += [xx]
    return result






def km_gwt_decode_contraction(term):                                        # \eta(p, q) = \delta(p, q) - \gamma(p, q)
    return km_gwt_decode_delta_s(km_gwt_decode_eta_s([term]))



def km_gwt_decode_eta_s(terms):              
    result = []
    for xx in terms:
        result += km_gwt_decode_eta(xx)
    return result


def km_gwt_decode_eta(term):
    counter = 0
    for xx in term.coefficient.matElement:
        if(xx.name == 'eta'):
            counter += 1
    if(counter == 0):
        return [term]
    else:
        result = [copy.deepcopy(term)]                                         # counter is # of 'eta' in term, so it should be the same for all derivative terms
        for i in range(counter):
            result = km_gwt_decode_eta_1_s(result)                             # km_gwt_so_decode_eta_1_s() replace ONE 'eta' mat in every term by delta(p, q) - gamma(p, q)
        return result


def km_gwt_decode_eta_1_s(terms):
    result = []
    for xx in terms:
        result += km_gwt_decode_eta_1(xx)
    return result


def km_gwt_decode_eta_1(term):
    flag = 0
    others = []
    flag2 = 1
    for xx in term.coefficient.matElement:
        if(xx.name == 'eta' and flag2 == 1):                                   # since we only replace ONE 'eta', we use flag2 to monitor
            flag += 1
            eta = xx
            flag2 = 0
        else:
            others += [xx]
    if(flag == 0):
        printterm(term)
        pdb.set_trace()
        return []
    elif(system_switch['use spin orbital'] == 1):                   # eta(p, q) = delta(p, q) - lambda(p, q)
        if(eta.matUpperIndicees[0].type in ['p', 'P', 'e', 'E', 'a'] and eta.matLowerIndicees[0].type in ['h', 'H', 'm', 'M', 'i']):
            if(system_switch['EXCLUDE Particle Ind from CU'] == 0):                                  # in this case, 'delta' == 0
                return [Term(Coefficient(-1.0*term.coefficient.const, others + [MatElement('lambda1', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators)]
            else:
                return []
        elif(eta.matUpperIndicees[0].type in ['h', 'H', 'm', 'M', 'i'] and eta.matLowerIndicees[0].type in ['p', 'P', 'e', 'E', 'a']):
            if(system_switch['EXCLUDE Particle Ind from CU'] == 0):                                  # in this case, 'delta' == 0
                return [Term(Coefficient(-1.0*term.coefficient.const, others + [MatElement('lambda1', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators)]
            else:
                return []
        elif(system_switch['EXCLUDE Particle Ind from CU'] == 1 and (eta.matUpperIndicees[0].type in ['p', 'P', 'e', 'E', 'a'] or eta.matLowerIndicees[0] in ['p', 'P', 'e', 'E', 'a']) ):
            return [Term(Coefficient(term.coefficient.const, others + [MatElement('delta', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators)]
        else:
            return [Term(Coefficient(term.coefficient.const, others + [MatElement('delta', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators),
                    Term(Coefficient(-1.0*term.coefficient.const, others + [MatElement('lambda1', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators)]
    elif(system_switch['use spin orbital'] == 0):                   # eta(p, q) = 2 * delta(p, q) - Gamma(p, q)
        if(eta.matUpperIndicees[0].type in ['p', 'P', 'e', 'E', 'a'] and eta.matLowerIndicees[0].type in ['h', 'H', 'm', 'M', 'i']):
            if(system_switch['EXCLUDE Particle Ind from CU'] == 0):                                  
                return [Term(Coefficient(-1.0*term.coefficient.const, others + [MatElement('Lambda1', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators)]
            else:
                return []
        elif(eta.matUpperIndicees[0].type in ['h', 'H', 'm', 'M', 'i'] and eta.matLowerIndicees[0].type in ['p', 'P', 'e', 'E', 'a']):
            if(system_switch['EXCLUDE Particle Ind from CU'] == 0):          
                return [Term(Coefficient(-1.0*term.coefficient.const, others + [MatElement('Lambda1', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators)]
            else:
                return []
        elif(system_switch['EXCLUDE Particle Ind from CU'] == 1 and (eta.matUpperIndicees[0].type in ['p', 'P', 'e', 'E', 'a'] or eta.matLowerIndicees[0] in ['p', 'P', 'e', 'E', 'a']) ):
            return [Term(Coefficient(term.coefficient.const*2, others + [MatElement('delta', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators)]
        else:
            return [Term(Coefficient(term.coefficient.const*2, others + [MatElement('delta', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators),
                    Term(Coefficient(-1.0*term.coefficient.const, others + [MatElement('Lambda1', eta.matUpperIndicees, eta.matLowerIndicees)]), term.uOperators)]
    else:
        set_trace()
        return


def km_gwt_decode_delta_s(terms):
    result = []
    for xx in terms:
        result += km_gwt_decode_delta(xx)
    return result


def km_gwt_decode_delta(term):                                                 # decoding delta will transform 'delta'
    counter = 0
    for xx in term.coefficient.matElement:
        if(xx.name == 'delta'):
            counter += 1
    if(counter == 0):
        return [term]
    else:
        result = [copy.deepcopy(term)]                                         # counter is # of delta mats, which should be the same for all derivative terms
        for i in range(counter):
            result = km_gwt_decode_delta_1_s(result)                           # mn_gwt_decode_delta_1_s() replace 'delta' mat in every term ONLY ONCE!
        result = changematname_s(result, "XXX", "delta")
        return result


def km_gwt_decode_delta_1_s(terms):
    result = []
    for xx in terms:
        result += km_gwt_decode_delta_1(xx)
    return result

def km_gwt_decode_delta_1(termm):                                               # replace ONE delta matrix by equating one index to another index
    flag = 0
    flag2 = 1
    others = []
    term = copy.deepcopy(termm)
    for xx in term.coefficient.matElement:
        if(xx.name == 'delta' and flag2 == 1):                                 # since we only replace ONE delta, we use flag2 to monitor
            flag += 1
            delta = xx
            flag2 = 0
        else:
            others += [xx]
    if(flag == 0):
        printterm(term)
        pdb.set_trace()
        return
    else:
        if(delta.matUpperIndicees[0].type in ['p', 'P', 'e', 'E', 'a'] and delta.matLowerIndicees[0].type in ['h', 'H', 'm', 'M', 'i']):
            return []                                                          # due to km_gwt, we might meet this case  
        elif(delta.matLowerIndicees[0].type in ['p', 'P', 'e', 'E', 'a'] and delta.matUpperIndicees[0].type in ['h', 'H', 'm', 'M', 'i']):
            return []
        elif(delta.matLowerIndicees[0].active == -1 and delta.matUpperIndicees[0].active == 1):
            return []
        elif(delta.matLowerIndicees[0].active == 1 and delta.matUpperIndicees[0].active == -1):
            return []
        elif(delta.matLowerIndicees[0].fix == 1 and delta.matUpperIndicees[0].fix == 1):   # deal with this special case
            delta.name = 'XXX'
            return [term]
        else:
            [ind1, ind2] = generalIndex(delta.matUpperIndicees[0], delta.matLowerIndicees[0])
            term = Term(Coefficient(term.coefficient.const, others), term.uOperators)
            return [reset(reset(term, [ind1], [ind2]), [ind2], [ind2])]        # when ind1's information (e.g., .att) is passed to ind2 while ind1 is not in 'term', without 
                                                                               # doing reset(.., ind2, ind2) the informaton is lost  







#*******************************************************************************
#  replace 2-dody cumulants by 1-body cumulants and 2-body RDM or the opposite *
#*******************************************************************************
            
    





def replace_2b_cumu_with_dm_s(terms):              

    """ for spin orbital: replace lambda(p,q; r,s) by \gamma(p,q; r,s) - \lambda(p;r) \lambda(q;s) + \lambda(p;s)\lambda(q,r);
        for spatial orbital: replace Lambda(p,q; r,s) by  \Gamma(p,q; r,s) - \Lambda(p;r) \Lambda(q;s) + 0.5 \Lambda(p;s)\Lambda(q,r);
        the goal is to see if number of terms can be reduced,    """

    result = []
    for xx in terms: result.extend(replace_2b_cumu_with_dm(xx))
    return result



def replace_2b_cumu_with_dm(term):
    def replace_2b_cumu_with_dm_1_s(terms):
        result = []
        for xx in terms: result.extend(xx.replace_2b_cumu_with_dm_1())
        return result
    if system_switch['use spin orbital'] == 1: CumulantName = 'lambda'
    else: CumulantName = 'Lambda'
    counter = 0
    for xx in term.coefficient.matElement:
        if(xx.name == CumulantName and len(xx.matUpperIndicees) == 2):
            counter += 1
    if(counter == 0):
        return [term]
    else:
        result = [copy.deepcopy(term)]                     
        for i in range(counter):
            result = replace_2b_cumu_with_dm_1_s(result)        
        return result










def replace_2b_dm_with_cumu_s(terms):              

    """ for spin orbital: replace gamma(p,q; r,s) by \lambda(p,q; r,s) + \lambda(p;r) \lambda(q;s) - \lambda(p;s)\lambda(q,r);
        for spatial orbital: replace Gamma(p,q; r,s) by  \Lamba(p,q; r,s) + \Lambda(p;r) \Lambda(q;s) - 0.5 \Lambda(p;s)\Lambda(q,r);
        the goal is to see if number of terms can be reduced,    """

    result = []
    for xx in terms: result.extend(replace_2b_dm_with_cumu(xx))
    return result



def replace_2b_dm_with_cumu(term):
    def replace_2b_dm_with_cumu_1_s(terms):
        result = []
        for xx in terms: result.extend(xx.replace_2b_dm_with_cumu_1())
        return result
    if(system_switch['use spin orbital'] == 1): DMName = 'gamma'
    else: DMName = 'Gamma'
    counter = 0
    for xx in term.coefficient.matElement:
        if(xx.name == DMName and len(xx.matUpperIndicees) == 2): counter += 1
    if(counter == 0): return [term]
    else:
        result = [copy.deepcopy(term)]                     
        for i in range(counter): result = replace_2b_dm_with_cumu_1_s(result)        
        return result






def replace_12b_sf_dm_with_delta_s(terms):              

    """ For spatial orbital: replace Gamma(p,q; r,s) by  4 delta(p,r) * delta(q, s) -  2 delta(p, s) delta(q, r)
        Gamma(p, q) = 2 delta(p, q); this function works regardless of # of Gammas  """
    result = []
    for xx in terms: result.extend(replace_12b_sf_dm_with_delta(xx))
    return result



def replace_12b_sf_dm_with_delta(term):
    if system_switch['use spin orbital'] != 0: set_trace()
    def replace_12b_sf_dm_with_delta_1_s(terms):
        result = []
        for xx in terms: result.extend(xx.replace_12b_sf_dm_with_delta_1())
        return result
    DMNames = ['Gamma', 'Lambda1']
    counter = 0
    for xx in term.coefficient.matElement:
        if xx.name in DMNames and len(xx.matUpperIndicees) in [1, 2]: counter += 1
    if(counter == 0): return [term]
    else:
        result = [copy.deepcopy(term)]                     
        for i in range(counter): result = replace_12b_sf_dm_with_delta_1_s(result)        
        return result









def replace_2b_cu_with_specified_inds_term(ta, specify_inds):       
    """
        The current (new) version works for both spin-free and spin-orbital.
        replace those 2-body-lambda, which contain the 'specify_inds', by 2-body-DM and 1-body-DM whose upperindices are purely from one term and lowerindicees are purely from the other term.
        the reason is that the 2-body-lambda terms with those binary 1-1 contractions will combine to form a 2-body-DM; therefore the number of terms is
        reduced. Looking in a different perspective, a factorziation is achieved. 
        After calling solvetwoop(), we immediately discard terms of quandratic 2-b-cumulants, because discartion needs to use cumulants directly; therefore, if
        the replacement occurs, some cumulants are converted to gamma, then we can not easily discard those terms using a simple criterion. 
    """
    if 1:
        if system_switch['use spin orbital'] == 1: CumuName = 'lambda'
        else: CumuName = 'Lambda'
        if 1:
            if len(specify_inds):
                valid = 0
                for yy in ta.coefficient.matElement:
                    if(yy.name == CumuName and len(yy.matUpperIndicees) == 2):
                        if(len(OVERLAP(specify_inds, yy.matUpperIndicees + yy.matLowerIndicees)) == len(specify_inds)): 
                            valid = 1
                            break
            else: valid = 1
            if valid == 0: return [ta]
            else: return replace_2b_cumu_with_dm(ta)
          

    
def replace_2b_cu_with_specified_inds_terms(Ts, specify_inds):  
    result = []
    for xx in Ts: result.extend(replace_2b_cu_with_specified_inds_term(xx, specify_inds))
    return result















#--------------------------------------------------------------------------------------------------------------------------------------

#                                                  Tools Library

#--------------------------------------------------------------------------------------------------------------------------------------






class ToolsSeparator:
    def __init__(self, start = "Y"):
        self.start = type










TOOLS_LIBRARY = { 'Avoid Index Overlapping':  ['inc_index_avoid_overlap',
                                               'max_ind_num',
                                               'max_ind_num_term'],
                           'density_matrix':  ['mc_density',
                                               'non_0_density'],
                      'Coding and Decoding':  ['code_term',
                                               'code_term_s',
                                               'right_config'],
                               'Miscellous':  ['read_function',
                                               'find_first_available_ind_num',
                                               'filter',
                                               'listtype',
                                               'factorial',
                                               'removezero',
                                               'dress_a',
                                               'opposite_spin_s',
                                               'opposite_spin',
                                               'dress',
                                               'rev_term',
                                               'rev_terms',
                                               'dress_rev',
                                               'normalize',
                                               'convert_a_type',
                                               'pickco',
                                               'varyco',
                                               'changeco',
                                               'con_max',
                                               'give_offset_s',
                                               'delete',
                                               'deletes',
                                               'zerocoeff',
                                               'delete_zerocoeff',
                                               'interchange',
                                               'conjugate',
                                               'conjugates',
                                               'makeIndexSet',
                                               'combineindlist',
                                               'degeneracy',
                                               'purify',
                                               'keep_1',
                                               'give_full_type',
                                               'give_full_type_op',
                                               'give_full_type_mats',
                                               'find_contracted_inds',
                                               'give_indlist',
                                               'indlist_equal',
                                               'check_repeat_mat',
                                               'check_repeat_mats',
                                               'check_repeat_term',
                                               'check_repeat_terms',
                                               'density_times_2',
                                               'perturbation_order',
                                               'filter_perturbation_order',
                                               'select_perturbation_order',
                                               'makeup_external_s',
                                               'give_all_unique_ind',       # give all unique indices of a term
                                               'give_all_ind',              # give all indices of a term
                                               'combine_2_mat_to_1',
                                               'so_triple_coeff_s',
                                               'triple_coeff_s',
                                               'Four_DM_replace',
                                               'replace_one_mat_term_s',
                                               'keep_h_connectivity'],
                                'partition':  ['partition_int',
                                               'group_partition_int'],
                                   'search':  ['findposition',
                                               'findfree',
                                               'summationindices',
                                               'findexternal'],
                            'Judge: 0 or 1':  ['overlap',
                                               'OVERLAP',
                                               'two_terms_connect',
                                               'both_in_one_op'],
                                   'parity':  ['transpairty'],
                                    'count':  ['numtype',
                                               'appear_time',
                                               'num_appear_time',
                                               'countmt',
                                               'counttype',
                                               'counttype_mat',
                                               'counttype_op',
                                               'counttype_mats',
                                               'countmt_2_va'],
                          'permu and combi':  ['antisymmetrize_pq',
                                               'antisymmetrize_pq_s',
                                               'permute12',
                                               'permute12_s',
                                               'combi',
                                               'rmEl',
                                               'permu',
                                               'permu_general',
                                               'combi_lists',
                                               'permu_2',
                                               'permu_1_va',
                                               'unique_section_divide', 
                                               'permu_unique_section_divide',
                                               'inter_out_permu_unique_section_divide',
                                               'inter_permu_list',
                                               'combi_pick_one'],
                          'symmetry_number':  ['symmetry_number',
                                               'symmetry_number_2',
                                               'symmetry_number_3',
                                               'coeff_from_mats'],
                                 'specials':  ['setchar',
                                               'setchars',
                                               'pickupchar',
                                               'set',
                                               'setterms',
                                               'setfix',
                                               'setfixs',
                                               'reset',
                                               'resets',
                                               'inc_index',
                                               'inc_index_s',
                                               'break_g',
                                               'expand_general_indices',
                                               'set_inactive',
                                               'set_inactive_s',
                                               'set_term',
                                               'set_terms',
                                               'D_2_OP',
                                               'D_2_OP_1',
                                               'delete_x_s',
                                               'delete_x'],
                  'product and commutator':  ['formproduct',
                                               'product',
                                               'normal_product',
                                               'contraction_term_term',
                                               'steom_contraction_term_term',
                                               'contraction_term_terms',
                                               'steom_contraction_term_terms',
                                               'contraction_terms_term',
                                               'steom_terms_term',
                                               'contraction_terms_terms',
                                               'steom_contraction_terms_terms',
                                               'commutator',
                                               'able_full_contraction',
                                               'full_contraction_commutator'],
                                'selecting':  ['discard_3d_term_s',
                                               'discard_terms_have_mat',
                                               'select_terms_rank_n',
                                               'pickup',
                                               'final_pick()',
                                               'pickrank()',
                                               'choosetype()',
                                               'give_type()',
                                               'give_mattype()',
                                               'pick_type()',
                                               'choose_mat_name()',
                                               'choose_terms_matname()',
                                               'choose_terms_matnames()',
                                               'discard_terms_matname()',
                                               'choose_preference()',
                                               'pickup_terms_of_n_mats()',
                                               'select_ex_order()', 
                                               'independent_mat',
                                               'select_up_2_n_B_D()',
                                               'remove_core_cumu()'],
                  'comparing & reordering ':  ['fixorder()',
                                               'fixorder_s()',
                                               'term_prefer_mats()',
                                               'terms_prefer_mats()',
                                               'finalorder()',
                                               'dif()',
                                               'indlisteq()',
                                               'gvhighestpos()',
                                               'gvSmallest()',
                                               'give_highest_indnum()',
                                               'rank()',
                                               'rank_indices_by_num()',
                                               'pick_1_name()',
                                               'reorder_inds_hmpe()',
                                               'reorder()',
                                               'reorder_s()',
                                               'max_match()'],
                      'term(s) <--> string':  ['parseEndPoint()',
                                               'term2string()',
                                               'string2term()'],
                            'texform title':  ['title_1()',
                                               'build_title()',
                                               'build_title_3()',
                                               'build_title_2()'],
                               'changename':  ['changematname()',
                                               'changematname_s()',
                                               'changename()',
                                               'changename_s()',
                                               'changeindname()',
                                               'changeindname_s()']}





#*************************************************************************************
#  Avoid Index Overlapping (introduction of 'DEN' and 'Cbar' necessitates this step) *
#*************************************************************************************





def inc_index_avoid_overlap(terms1, terms2):             
    """ compare with terms2, increase the num's of indices in terms1 to avoid index overlap """

    P = max_ind_num(terms2, 'p')
    H = max_ind_num(terms2, 'h')
    M = max_ind_num(terms2, 'm')
    E = max_ind_num(terms2, 'e')
    G = max_ind_num(terms2, 'g')
    terms1 = inc_index_s(terms1, 'p', P)
    terms1 = inc_index_s(terms1, 'e', E)
    terms1 = inc_index_s(terms1, 'h', H)
    terms1 = inc_index_s(terms1, 'm', M)
    terms1 = inc_index_s(terms1, 'g', G)
    return terms1



def max_ind_num(terms, type):          
    """ return the maximum 'num' of indices of type 'type' in terms """

    result = 0
    for xx in terms:
        nu = max_ind_num_term(xx, type)
        if nu > result: result = nu
    return result


def max_ind_num_term(term, type):      
    """ return the maximum 'num' of indices of type 'type' in term """
    inds = give_all_ind(term)
    result = 0
    for xx in inds:
        if(xx.type == type):
            if(int(xx.num) > result): result = int(xx.num)
    return result




                                                     

#**************************************************
#                density-matrix                   *
#**************************************************





def non_0_density(terms):          

    """ if any index of operator is non hole, then the corresponding density matrix element will be zero, we eleminate those terms """

    result = []                       
    for xx in terms:
        if(len(xx.uOperators[0].upperIndicees) == 0):
            result += [xx]
        else:
            ins = xx.uOperators[0].upperIndicees + xx.uOperators[0].lowerIndicees
            flag = 1
            for yy in ins:
                if(yy.type in ['p', 'P', 'e', 'E']):                           # this is better than keeping 'm' and 'h', because sometines 'g' indices are there
                   flag = 0
                   break
            if(flag == 1):
                result += [xx]
    return result



def mc_density(terms1):           

    """ transform operator to density matrices, used in multireference cases """

    result = []
    if(system_switch['normal order'] in ['T', 'V']):
        terms = non_0_density(terms1)      # Pay attenstion to the order of upp and low !!!!!!
    elif(system_switch['normal order'] in ['KM']):
        terms = filter(terms1, 0)
    else:
        print "\n normal order undefined \n"
        pdb.set_trace()
    for xx in terms:
        if(xx.uOperators[0].upperIndicees == []):
            result += [xx]
        elif(system_switch['normal order'] in ['V', 'T']):
            namedir= {1: 'd1', 2: 'd2', 3: 'd3'}
            xx.coefficient.matElement += [MatElement(namedir[len(xx.uOperators[0].upperIndicees)], xx.uOperators[0].upperIndicees, xx.uOperators[0].lowerIndicees)]
            xx.uOperators[0] = Uoperator([], [])
            result += [xx]
        elif(system_switch['normal order'] == 'KM'):
            pass
        else:
            print "\n normal order undefined in mc_density()\n"
            pdb.set_trace() 
    if(system_switch['DM coefficient'] == 'normal'):
        pass
    elif(system_switch['DM coefficient'] == '1/n!' and (system_switch['Max DM rank'] in [2,3]) ):
        result = density_times_2(result)
#   elif(system_switch['DM coefficient'] == '1/n!' and (system_switch['Max DM rank'] in [4]) ):
#       print "\n treat Four_DM differently \n"
#       pass
    else:
        print "\n\n                    Unspecified !"
        pdb.set_trace()
    return result






#**************************************************
#               Coding and Decoding               *
#**************************************************





def code_term(term):         

    """ we code the information of position for external indices (or indices on operator) using Index.att """

    NATURE = term.uOperators[0]
    if(len(NATURE.upperIndicees) == 0):
        pass
    elif(len(NATURE.upperIndicees) == 1):
        NATURE.upperIndicees[0].att = 'U1'
        NATURE.lowerIndicees[0].att = 'L1'
    elif(len(NATURE.upperIndicees) == 2):
        NATURE.upperIndicees[0].att = 'U1'
        NATURE.upperIndicees[1].att = 'U2'
        NATURE.lowerIndicees[0].att = 'L1'
        NATURE.lowerIndicees[1].att = 'L2'
    else:
        print "\n\n\n three body operator? we are in code_term()"
        pdb.set_trace()
    inds = NATURE.upperIndicees + NATURE.lowerIndicees
    counter = 0
    for yy in term.coefficient.matElement:
        for zz in yy.matUpperIndicees + yy.matLowerIndicees:
            for zzz in inds:
                if(zz.equal(zzz)):
                    zz.att = zzz.att
                    counter += 1
                    break
    if(counter != len(inds)):        # for safety, we do a simple check       
        print "\n\n\n problem2 in code_term!\n\n"
        pdb.set_trace() 
    else:
        return term     


    

def code_term_s(terms):
    res = []
    for xx in terms:
        res += [code_term(xx)]
    return res


def right_config(uop):                 

    """ used accompanying multilevel condensing, to adjust the the indicees in operators to right positions, (these indices is coded in some way) """
    """ currently I assume the Uoperaor can't be null """

    R = Uoperator()    
    dic_right_config = {}
    for xx in uop.upperIndicees + uop.lowerIndicees:
        dic_right_config[xx.att] = xx
    if(len(dic_right_config.keys()) == 2):
        R.upperIndicees = [dic_right_config['U1']]   
        R.lowerIndicees = [dic_right_config['L1']]
    elif(len(dic_right_config.keys()) == 4):
        R.upperIndicees = [dic_right_config['U1'], dic_right_config['U2']] 
        R.lowerIndicees = [dic_right_config['L1'], dic_right_config['L2']]
    else:
        print "\n\n problem1 in right_config()!\n"
        pdb.set_trace()
    return R





#**************************************************
#                  Miscellous                     *
#**************************************************




def find_first_available_ind_num(term, indtype):                      

    """ find the first available/unused number of ind of type 'indtype', return an integer. e.g. term = V^i1_a1 * E^a1_i1, indtype = 'a', return 2 """  

    indlist = give_indlist(term)                                     
    counts = []
    for xx in indlist:
        if(xx.type == indtype): 
            counts += [int(xx.num)]
    flag = 0
    for count in range(1000):
        countt = count + 1
        if(countt not in counts):
            flag = 1
            break
    if(flag == 1):
        return countt
    else:
        print "\n\n Not found. Please check find_first_available_ind_num()\n"
        pdb.set_trace()
        return



def find_20_available_ind_num(term, indtype): 
    indlist = give_indlist(term)             
    counts = []
    for xx in indlist:
        if(xx.type == indtype):
            counts += [int(xx.num)]
    flag = 0
    for count in range(1000):
        countt = count + 1
        if(countt not in counts):
            counts += [countt]
            flag += 1
            if(flag >= 20):
                break
    if(flag == 20):
        return countt
    else:
        print "\n\n Not found. Please check find_first_available_ind_num()\n"
        pdb.set_trace()
        return




def filter(terms, n):            
    """ filter terms of rank over n """
    result = []
    for a in terms:  
        if(len(a.uOperators) == 0): 
            printterm(a)
            print "\n unexpected; normally there should be at least one operator without indices \n"
            set_trace()
        if len(a.uOperators[0].upperIndicees) <= n: result.append(a) 
    return result


def listtype(inds):                   
    """ return the type of list of ind """
    result = ""
    for x in inds: result += x.type
    return result
    



def factorial(n):
    res = 1.0
    i = 1
    if n >= 0 and n < 1: return res
    else:
        while(i < n + 1):
            res = res*i
            i = i + 1
        return res




def removezero(A):                     # remove constants in terms, employed in commutator()
    res = []
    for x in A:
        if len(x.uOperators[0].upperIndicees) == 0: pass
        else: res = res + [x]
    return res



def dress_a(term):
    return dress([term])[0]



def opposite_spin_s(terms):
    """ convert the spin of every index to the opposite spin """
    result = []
    for xx in terms:
        result.append(opposite_spin(xx))
    return result



def opposite_spin(term):
    dic = {'p': 'P', 'P': 'p',
           'h': 'H', 'H': 'h',
           'm': 'M', 'M': 'm',
           'e': 'E', 'E': 'e'}
    term2 = copy.deepcopy(term)
    for xx in term2.coefficient.matElement:
        for yy in xx.matLowerIndicees + xx.matUpperIndicees:
            yy.type = dic[yy.type]
    return term2    

              
    
def dress(terms):                      
    """ renumber indices of terms, first operator then coefficient """
    res = []
    counter = 0
    for y in terms:
        x = copy.deepcopy(y)
        ind = makeIndTerm(x)
        ind = ind.renumber()
        x = arrangeIndicees(ind, x)
        res.append(x)
        counter += 1
    return res



def pickco(terms):                     # collect the coefficients of terms, i.e., set the operator to be empty
    res = []
    for xx in terms:
        yy = copy.deepcopy(xx)
        yy.uOperators = [Uoperator([], [])]
        res.append(yy) 
    return res



def multiply_co(co, terms):                 # multiply each term by a constant "co"
    result = []
    for xx in terms:
        yy = copy.deepcopy(xx)
        yy.coefficient.const *= co
        result.append(yy) 
    return result



def divide_co(co, terms):                 # divide each term by a constant "co"
    result = []
    for xx in terms:
        yy = copy.deepcopy(xx)
        yy.coefficient.const /= co
        result.append(yy) 
    return result



def changeco(co, terms):                
    result = []
    for xx in terms:
        yy = copy.deepcopy(xx)
        yy.coefficient.const = co
        result.append(yy) 
    return result




def give_offset_s(nums, listt):
    result = []
    for xx in nums:
        result += [listt.index(xx)]
    return result




def delete(list, offset):              # delete the elment in list of offset "offset", if it exceeds the length of list, print "swap"
    if(offset + 1 > len(list)):
        print "swap!!"
        return
    else:
        if(offset < len(list) - 1): return list[:offset] + list[offset+1:]
        else: return list[:offset]




def deletes(listt, offsets):           # delete several elements, characterized offset, from a list
    result = []
    for n in range(len(listt)):
        if(n in offsets): pass
        else: result += [listt[n]]
    return result





def conjugate(term1):        
    """ produce the conjugate term, by interchanging upper and lower indices of operator, here we take the matrixelements to be symmetric, so nothing is done on them """
    term = copy.deepcopy(term1)       
    uppi = copy.deepcopy(term.uOperators[0].upperIndicees)
    lowi = copy.deepcopy(term.uOperators[0].lowerIndicees)
    term.uOperators[0].upperIndicees = lowi
    term.uOperators[0].lowerIndicees = uppi
    return term



def conjugates(terms):                 # produce the conjugate terms
    result = []
    for xx in terms:
        zz = conjugate(xx)
        result.append(zz) 
    return result


def invert_upp_low(term):              # interchange upper and lower indices of matelements of a term
    mats = []
    for xx in term.coefficient.matElement:
        mats += [MatElement(xx.name, xx.matLowerIndicees, xx.matUpperIndicees)]
    return Term(Coefficient(term.coefficient.const, mats), term.uOperators)



def invert_upp_low_s(terms):
    result = []
    for xx in terms:
        result += [invert_upp_low(xx)]
    return result 



def makeIndexSet(upperInd, lowerInd):
    """ makes an object of the class IndexSet from the labels """
    iSet = []
    if(len(upperInd) == 0): return IndexSet([])
    else:
        for b in range(len(upperInd)):
            iPair = IndexPair(upperInd[b], lowerInd[b], 1)
            iSet.append(iPair)
        return IndexSet(iSet)



def combineindlist(indlist1, indlist2):# find the combination of two lists of index, discard repeated ones
    L1 = indlist1
    for x in indlist2:
        if not x.inArray(L1): L1.append(x) 
    return L1




def degeneracy(type):                  # this function calculate the multiplying factor of one type of operator from type, i.e. calculate multiplicity
    if len(type) == 0: return 1
    else:
        L = []
        uppers = type[:len(type)/2]
        lowers = type[len(type)/2:]
        for i in range(len(uppers)): L += [uppers[i] + lowers[i]]
        LL = permu(len(L), L)
        j = 0
        for xx in LL:
            if xx == L: j += 1
        return j



def purify(list):                      # delete repeated element (keep only once)
    result = []
    for xx in list:
        flag = 0
        for x in result:
            if xx == x:
                flag = 1
                break
        if(flag == 0): result.append(xx) 
    return result



def give_full_type(term):              # given a term, return a list consiting of its index types, if some index is fixed, we also include its 'fix', 
    result = []                        # because fixed indices are special, this function will be used in canocomp_va() 
    for xx in term.coefficient.matElement:
        result += [xx.name]
        for yy in xx.matUpperIndicees + xx.matLowerIndicees:
            if(yy.fix == 1): result += [yy.type + yy.num] 
            else: result += [yy.type]
    for zz in term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees:
        if(zz.fix == 1): result += [zz.type + yy.num]
        else: result += [zz.type]
    return result



def give_full_type_op(op):
    result = []
    for zz in op.upperIndicees + op.lowerIndicees:
            result += [zz.type]
    return result


def give_full_type_mats(mats):
    result = []
    for xx in mats:
        for yy in xx.matUpperIndicees + xx.matLowerIndicees: result += [yy.type]
    return result


def find_contracted_inds(inds):        # find contracted indices in 'inds'
    res1 = []
    res2 = []
    for xx in inds:
        if(xx.inArray(res1)):
            if(xx.inArray(res2)):
                print "\n\n            Some index appearing three times.  find_contracted_inds@real.py\n\n"
                pdb.set_trace()
            else:
                res2 += [xx]
        else:
            res1 += [xx]
    return res2


def give_all_ind(term):                       # give all indices in term
    result = term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees
    for xx in term.coefficient.matElement:
        result += xx.matUpperIndicees + xx.matLowerIndicees
    return result




def give_indlist(term):                # give indices in a term, repeated ones only appearing once
    result = []
    for xx in term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees:
        if not xx.inArray(result): result.append(xx) 
    for yy in term.coefficient.matElement:
        for zz in yy.matUpperIndicees + yy.matLowerIndicees:
            if not zz.inArray(result): result.append(zz)
    return result 



def uniqueinds_from_indlist(inds):    
    """ return all unique indices """ 
    result = []
    for zz in inds:
        if(zz.inArray(result) == 0): result.append(zz) 
    return result 



def pick_type_ind_from_inds(inds, ttype):
    """ return those inds of type 'type' """
    res = []
    for xx in inds:
        if xx.type == ttype: res.append(xx)
    return res




def indlist_equal(indL1, indL2):
    result = 1
    if(len(indL1) != len(indL2)): result = 0
    else:
        for xx in range(len(indL1)):
            if(indL1[xx].equal(indL2[xx]) == 0): result = 0
    return result



def check_repeat_mat(mat):             # return 1 if finding repeated indices in a matrix element
    res = 0
    JJ = []
    for xx in mat.matUpperIndicees + mat.matLowerIndicees: JJ += [xx.type + xx.num]
    for i in range(len(JJ)-1):
        for j in  range(len(JJ)-1-i):
            if(JJ[i] == JJ[i + j +1]):
                res = 1                # find repeated indices
                printMatElement(mat)
                set_trace()
                break
    return res



def check_repeat_mats(mats):
    res = 0
    for xx in mats:
        if(check_repeat_mat(xx)): 
            res = 1
            break
    return res



def check_repeat_term(term):
    Y = term.coefficient.matElement
    for xx in term.uOperators:
        Y += [MatElement('YYYY', xx.upperIndicees, xx.lowerIndicees)]
    return check_repeat_mats(Y)


    
def check_repeat_terms(terms):
    for xx in terms:
        if(check_repeat_term(xx)):
            pdb.set_trace()
            printterm(xx)
            pdb.set_trace()
    return                   



def density_times_2(terms):
    result = []
    for xx in terms:
        for yy in xx.coefficient.matElement:
            if(yy.name in['nd', 'd1', 'd2', 'd3', 'd4', 'x1', 'x2', 'x3']):
                if(len(yy.matUpperIndicees) == 1): pass
                elif(len(yy.matUpperIndicees) == 2):
                    xx.coefficient.const *= 2.0
                    break
                elif(len(yy.matUpperIndicees) == 3):    
                    xx.coefficient.const *= 6.0                                     
                    break
                elif(len(yy.matUpperIndicees) == 4):
                    xx.coefficient.const *= 24.0
                    break
                else:
                    print "\n condition undefined\n"
                    pdb.set_trace()
        result.append(xx) 
    return result




def makeup_external_s(terms, indices, atts):                                   # by comparing 'att', we identify the corresponding indices
    result = []                                                                # then by adding delta tensor to force the external indices
    for xx in terms: result.append(makeup_external(xx, indices, atts))         # to be 'indices'
    return result



def makeup_external(term, indices, atts):
    if(len(indices) != len(atts)):
        print "\n\n                    Error in makeup_external@real.py\n\n"
        print len(indices), len(atts)
        printterm(term)
        pdb.set_trace()
        return
    elif(len(indices) == 0): return copy.deepcopy(term)
    else:
        result = copy.deepcopy(term)
        for n in range(len(indices)): result = makeup_external_1_ind(result, indices[n], atts[n])
        return result



def makeup_external_1_ind(term, ind, attt):
    result = copy.deepcopy(term)
    exs = findexternalinds(result)
    flag = 0
    for xx in exs:
        if(xx.att == attt):
            flag += 1
            indd = xx
    if(flag == 0):
        print "\n\n                    Error: attribute lost somewhere\n\n"
        printterm(term)
        print attt
        for yy in exs:
            print yy.gvIndex()
            print yy.att
        pdb.set_trace()
    elif(flag > 1):
        print "\n\n                    Error: more than one index have the same attribute\n\n"
        printterm(term)
        print attt
        for yy in exs:
            print yy.gvIndex()
            print yy.att
        pdb.set_trace()
    else:
        if(indd.gvIndex() == ind.gvIndex()): return result
        else:
            result.coefficient.matElement += [MatElement('delta', [indd], [ind])]
            print attt
            print ind.gvIndex()
            printterm(term)
            printterm(result)
            pdb.set_trace()
            return result






def so_triple_coeff(term):                   # if a term contains a 't3' term, multiply the coefficient by 1/36
    for xx in term.coefficient.matElement:
        if(xx.name == 't3'): term.coefficient.const /= 36.0
    return term



def so_triple_coeff_s(terms):
    result = []
    for xx in terms: result.append(so_triple_coeff(xx))
    return result



def triple_coeff(termm):                      # if a term contains a 't3' term, multiply the coefficient by 1/6
    terma = copy.deepcopy(termm)
    for xx in terma.coefficient.matElement:
        if(xx.name == 't3'): terma.coefficient.const /= 6.0
    return terma 



def triple_coeff_s(terms):
    result = []
    for xx in terms: result += [triple_coeff(xx)]
    return result



def keep_h_connectivity(terms):        # only keep terms which are connceted to H     
    res = []
    for xx in terms:
        if(h_connectivity(xx) == 1):
            res += [xx]
    return res



def h_connectivity(termm):
    conn = 1
    find_h = 0
    debug = 0
    for xx in termm.coefficient.matElement:
        if(xx.name in ['f', 'v', 'nv', 'nf', 'kmf', 'kmv']):
            find_h += 1
            v_inds = xx.matUpperIndicees + xx.matLowerIndicees
    if(find_h != 1):
        print "\n  \n"
        printterm(termm)
        pdb.set_trace()
    else:
        if(len(termm.coefficient.matElement) > 1):
            for xx in termm.coefficient.matElement:
                if(xx.name in ['f', 'v', 'nv', 'nf', 'kmf', 'kmv']):
                    pass
                elif(xx.name in ['DEN', 'Cbar', 'delta', 'eta', 'Lambda1', 'lambda1', 'Lambda', 'lambda'] and (overlap(xx.matUpperIndicees + xx.matLowerIndicees, v_inds) == 1)):
                    v_inds += xx.matUpperIndicees + xx.matLowerIndicees
            for xx in termm.coefficient.matElement:
                if(xx.name in ['f', 'v', 'nv', 'nf', 'kmf', 'kmv', 'DEN', 'delta', 'Cbar', 'eta', 'Lambda1', 'lambda1', 'Lambda', 'lambda']): pass
                elif(overlap(xx.matUpperIndicees + xx.matLowerIndicees, v_inds) == 1): pass
                else:
                    conn = 0
                    break
        return conn






#**************************************************
#                   partition                     *
#**************************************************


def partition_int(M):

    """ partition an integer to all possible parts, e.g. 1-> [ [1] ], 2-> [ [2], [1,1] ], 3-> [ [3], [1,2], [2,1], [1,1,1] ] """

    result = []
    if(M == 0):
        print " in partition_int : M == 0"
        pdb.set_trace()
    elif(M == 1):
        result = [ [1] ]
    else:
        for ma in range(M):
            mb = ma + 1
            if(mb == M):
                result += [ [M] ]
            else:
                mc = M - mb
                part = partition_int(mc)
                for xx in part:
                    xx = [mb] + xx
                    result += [ xx ]
    return result


def group_partition_int(M):

    """ partition an integer to all possible group parts based in increasing section size, e.g. 1-> [ [1] ], 2 -> [ [2], [1,1] ], 3-> [ [3], [1,2], [1,1,1] ]"""
    """ [2,1] will not be included. To fold info compactly, we then just record the infor the number appearing in every non-empty section, again in """
    """ increasing section size. Then we code info in Uoperator, e.g., M = 3, we return [ op(3;1), op(1,2; 1,1), op(1; 3) ] """

    debug = 0
    if(debug):
        print "\n debug group_partition_int \n"
        pdb.set_trace()
    result = []
    unique_parts = group_partition_int_help(M, 0)
    for xx in unique_parts:
        unique_nums = []
        upp_inds = []
        low_inds = []
        for yy in xx:
            if(yy not in unique_nums):
                if(len(unique_nums) > 0):
                    if(yy <= unique_nums[-1]):             # check increasing order
                        print "\n in group_partition_int, number not increasing? print out unique_nums and yy:\n"
                        print unique_nums
                        print yy
                        print "\n"
                        pdb.set_trace()
                unique_nums += [yy]
        for zz in unique_nums:
            upp_inds += [zz]
            low_inds += [num_appear_time(zz, xx)]
        result += [Uoperator(upp_inds, low_inds)]
    return result


def group_partition_int_help(M, minn):
    result = []
    if(M == 0):
        print " in partition_int : M == 0"
        pdb.set_trace()
    if(M < minn):
        return []
    elif(M == minn):
        result = [ [minn] ]
    else:
        for ma in range(M):
            mb = ma + 1
            if(mb >= minn):
                if(mb == M):
                    result += [ [M] ]
                else:
                    mc = M - mb
                    part = group_partition_int_help(mc, mb)
                    for xx in part:
                        xx = [mb] + xx
                        result += [ xx ]
    return result







#**************************************************
#                   search                        *
#**************************************************





def findposition(num, numlist):        # return the position of num in numlist, used in trasnspairty()
    pos = -1
    for x in range(len(numlist)):
        if(num == numlist[x]):
            pos = x
    return pos




def findfree(ind):                     # find free indices in a indlist, (in action)
    res = []
    for xx in ind:
        if(xx.fix == 0):
            res = res + [xx]
    return res



def findfix(inds):                     # find fix indices in a indlist, (in action)
    res = []
    for xx in inds:
        if(xx.fix == 1):
            res = res + [xx]
    return res




def summationindices(term):            # so far unsed
    LIST = []
    help = []
    for xx in term.coefficient.matElement:
        for yy in xx.matUpperIndicees + xx.matLowerIndicees:
            if(yy.inArray(help)):
                if(yy.inArray(LIST)): pass
                else: LIST += [yy]
            else: help += [yy]
    return LIST



def findexternal(mats):                # find uncontracted indices
    flag = 0
    indlist = []
    for xx in mats:
        indlist += xx.matUpperIndicees + xx.matLowerIndicees
    for zz in indlist:
        if(appear_time(zz, indlist) < 2):
            flag = 1
            break
    return flag



def findexternalinds(term):
    indlist = give_all_ind(term)
    result = []
    for zz in indlist:
        if(appear_time(zz, indlist) == 1):
            result += [zz]
    return result



def find_ex_lists(list1, list2):
    result = []
    L = list1 + list2
    for zz in L:
        if(appear_time(zz, L) == 1):
            result += [zz]
    return result

    



            




#**************************************************
#                 Judge: 0 or 1                   *
#**************************************************





def overlap(list1, list2):

    """ determine if there are overlap indices between list1 and list2 """

    if(len(list1) == 0): return 0
    elif(len(list2) == 0): return 0
    else:
        result = 0
        for xx in list1:
            if(xx.inArray(list2)):
                result = 1
                break
        return result



def OVERLAP(list1, list2):

    """ give the overlap indices """

    if(len(list1) == 0): 
        return []
    elif(len(list2) == 0): 
        return []
    else:
        result = []
        for xx in list1:
            if(xx.inArray(list2)):
                result += [xx]
        return result



def two_terms_connect(term1, term2):   # whether the two terms are connected or not
    inds1 = give_all_ind(term1)
    inds2 = give_all_ind(term2)
    return overlap(inds1, inds2)



def both_in_one_op(A, B, LL):

    """ return 1 if (all numbers in A and B < LL OR all numbers in A and B >= LL); return 0 otherwise """

    A_from_op1 = []
    A_from_op2 = []
    B_from_op1 = []
    B_from_op2 = []
    for xx in A:
        if(xx < LL):
            A_from_op1 += [xx]
        else:
            A_from_op2 += [xx]
    if(len(A_from_op1) > 0 and len(A_from_op2) > 0):
        result = 0
    else:
        for xx in B:
            if(xx < LL):
                B_from_op1 += [xx]
            else:
                B_from_op2 += [xx]
        if(len(A_from_op1) == 0 and len(B_from_op1) == 0):
            result = 1
        elif(len(A_from_op2) == 0 and len(B_from_op2) == 0):
            result = 1
        else:
            result = 0
    return result





#**************************************************
#                   parity                        *
#**************************************************





def transpairty(listA, listB):         # determine the parity of transposition between two sequences, which is the positions of operators
    count = 0
    for i in range(len(listB)):
        position = findposition(listB[i], listA)
        if(position is not i):
            count = count + 1
            listA[position] = listA[i]
            listA[i] = listB[i]
    return pow(-1, count)





#******************************************
#                 count                   *
#******************************************





def numtype(term):                     # count # of various kinds of indices in a term, return # list
    ind = makeIndTerm(term)
    typelist = ["p", "e", "h", "m", "g"]
    list = []
    result = [0, 0, 0, 0, 0]
    for aa in ind.indexSets:
        for bb in aa.indexPairs:
            if(bb.upperIndex.inArray(list) == 0):
                result[typelist.index(bb.upperIndex.type)] += 1
                list = list + [bb.upperIndex]
            if(bb.lowerIndex.inArray(list) == 0):
                result[typelist.index(bb.lowerIndex.type)] += 1
                list = list + [bb.lowerIndex]
    return result




def appear_time(ind, list):            # to facilitate next function
    res = 0
    for xx in list:
        if(xx.equal(ind)):
            res += 1
    return res


def num_appear_time(num, list):
    res = 0
    for xx in list:
        if(xx == num):
            res += 1
    return res









def counttype(indlist):                # count # of each type of index in a indlist
    t1 = 0
    t2 = 0
    t3 = 0
    t4 = 0
    t5 = 0
    for x in indlist:
        if(x.type == "p"): t1 += 1
        elif(x.type == "e"): t2 += 1
        elif(x.type == "h"): t3 += 1
        elif(x.type == "m"): t4 += 1
        elif(x.type == 'g'): t5 += 1
        else:
            print "\n\n problem in counttype()  \n\n"
            pdb.set_trace()
    return [t1, t2, t3, t4, t5]




def counttype_mat(mat):                 
    """  # return a list, the first element is the list of number of each type of index in Upper; the second .... Lower ... """
    return counttype(mat.matUpperIndicees) + counttype(mat.matLowerIndicees)           


def counttype_op(term):
    return [counttype(term.uOperators[0].upperIndicees), counttype(term.uOperators[0].lowerIndicees)]



def counttype_mats(mats):
    result = []
    for x in mats:
        result += [counttype_mat(x)]
    return result






#**************************************************
#               permu and combi                   *
#**************************************************





def antisymmetrize_pq(termm, pq):

    """ antisymmetrize termm, return term(p,q) - term(q,p) """ 

    permu = reset(termm, pq, [ pq[1], pq[0] ])
    return [termm] + varyco(-1.0, [permu])
    


def antisymmetrize_pq_s(termms, pq):

    result = []
    for termm in termms:
        result += antisymmetrize_pq(termm, pq)
    return result


def antisymmetry_test(terms, pq):

    """ test whether 'terms' is antisymmetric """

    terms2 = resets(terms, pq, [pq[1], pq[0]])
    if(len(merge_all(terms + terms2)) == 0):
        return 1 
    else:
        return 0
       



def extract_unique_terms_from_explicit_antisym(terms, pq):

    """ suppose terms is explicitly antisymmetric, that is terms(p, q) can be written as X(p, q) - X(q, p).
        this function tries to find X so as to reduced the # of terms; this function should be used AFTER merge_all """
    
    result = []
    B = copy.deepcopy(terms)
    nn = len(terms)
    count = 0
    while (len(B) > 0):
        if(count == nn):                                                       # the maximum number of cycles is len(B)
            pdb.set_trace()
            return
        Y = [canonical_one_by_permu(reset(B[0], pq, [pq[1], pq[0]]))]
        if(len(merge_all(B[:1] + Y)) == 0):                                    # this means that X is itself antisymmetric
            result += varyco(0.5, B[:1])
            del B[0]
        else:
            findd = 0
            lenn = len(B)
            for xx in range(lenn)[1:]:
                if(len(merge_all([ B[xx] ] + Y)) == 0):
                    findd = 1
                    location = xx
                    break
            if(findd == 0):
                pdb.set_trace()                                                # not explicitly antisymmetrized?
                return
            else: 
                result += B[:1]
                del B[location] 
                del B[0]
        count += 1
    return result               





def permute12(term1, matname):         # assuming term has only one matelement which has name 'matname', and this matelment is a two-body one
    term = copy.deepcopy(term1)        # X^pq_rs -> X^pq_rs + X^qp_sr
    n = 0                              # In my case, matname is 'X'
    flag = 0
    for xx in term.coefficient.matElement:
        if(xx.name == matname):
            mat = xx
            n += 1
            flag = 1
    if(flag == 0):
        return [term1]
    else:
        if(n > 1):
            print "\n\n                P1@permute12.real!"
            pdb.set_trace()
        else:
            A = copy.deepcopy(mat.matUpperIndicees[0])
            B = copy.deepcopy(mat.matUpperIndicees[1])
            C = copy.deepcopy(mat.matLowerIndicees[0])
            D = copy.deepcopy(mat.matLowerIndicees[1])
            term = reset(term, mat.matUpperIndicees+mat.matLowerIndicees, [B, A, D, C])
            return [term1, term]



def permute12_s(terms, matname):
    res = []
    for xx in terms:
        res += permute12(xx, matname)
    return res



def combi(k, array):

        """ gives the combinations (ordered permutations) of any 'k' of objects from the array """

        result = []
        if(k == 1):
                for a in array:
                        result = result + [[a]]
                return result
        elif(k == len(array)):
                return [array]
        elif(k == 0):
                return []
        else:
                result = combi(k, array[1:])
                half = combi(k-1, array[1:])
                for b in half:
                        result = result + [[array[0]] + b]
                return result






def rmEl(a, arr):                      # used in permu(), unexamined

        """ removes the element at position 'a' from array 'arr' """

        if(len(arr) == 0): return []
        elif(arr.index(a) == len(arr) - 1):
            return arr[:arr.index(a)]
        else:
            return arr[:arr.index(a)] + arr[arr.index(a) + 1:]



def permu(k, array):

        """ gives the permutations of any 'k' of objects from the array """

        result = []
        if(k == 1):
                for a in array:
                        result = result + [[a]]
                return result
        elif(k == 0):
                return []
        else:
                for a in array:
                        half = permu(k-1,rmEl(a,array))
                        for b in half:
                                result = result + [[a] + b]
                return result




def permu_general(array):              # given an array, return all possible permuations
    k = len(array)
    result = []
    if k == 1:
        return [array]
    elif k == 0: return []
    else:
        for xx in range(len(array)):
            a = array[xx]
            half = permu(k-1,delete(array, xx))
            for b in half:
                result = result + [[a] + b]
        return result




def combi_lists(lists):              
    """ pick up one element from each A element and put them in a list, which is an element of B. 
        For example, if A = [[a1, a2], [b1, b2], [c1, c2]], then B = [ [a1, b1, c1], [a1, b1, c2], [a1, b2, c1], [a1, b2, c2], ...  ] """
    k = len(lists)    
    result = []
    if(k == 1):
         for a in lists[0]:
             result += [[a]]
         return result
    elif(k == 0):
         return []
    else:
        lista = lists[0]
        listss = delete(lists, 0)
        xyz = combi_lists(listss)
        for xx in lista:
            for yy in xyz:
                result += [[xx] + yy]
    return result                



def combi_lists_2(lists):              # obtaining all possible list by taking one element from each list
    k = len(lists)
    result = []
    if(k == 1):
         for a in lists[0]:            
             result += [a]
    elif(k == 0):
         return []
    else: 
        lista = lists[0]
        listss = delete(lists, 0)
        xyz = combi_lists_2(listss)
        for xx in lista:
            for yy in xyz:
                result += [xx + yy]
    return result



def permu_2(term1):                    # permute indexpairs in matelements and operator
    result = []
    term = copy.deepcopy(term1)
    indd = makeIndTerm(term)
    LL1 = []
    for xx in indd.indexSets:
        LL2 = []
        yy = xx.indexPairs
        LL3 = permu_general(yy)
        for zz in LL3:
            LL2 += [IndexSet(zz)]
        LL1 += [LL2]
    LL4 = combi_lists(LL1)
    for xxxx in LL4:
        result += [arrangeIndicees(Indicees(xxxx), term)]
    return result         




def permu_1_va(term1):                 # permute the same type of matelements
    print "\n This Fuction Should be Modified Before Being Used\n"
    pdb.set_trace()                    # I assume that this function is not in use, as there are certain types of matrix elements missing..
    result = []       
    term = copy.deepcopy(term1)
    ops = term.uOperators
    numm = term.coefficient.const
    LV = []
    LT1 = []
    LT2 = []
    LT3 = []
    LS1 = []
    LS2 = []
    LW = []
    LDEL = []
    LDEN = []
    LT = []
    for xx in term.coefficient.matElement:
        if(xx.name in ['f', 'v']):
            LV += [xx]
        elif(xx.name in ['d1', 'd2', 'd3', 'd4', 'x1', 'x2', 'x3']):
            LDEN += [xx]
        elif(xx.name == 'delta'):
            LDEL += [xx]
        elif(xx.name == 't1'):
            LT1 += [xx]
        elif(xx.name == 't2'):
            LT2 += [xx]
        elif(xx.name == 't3'):
            LT3 += [xx]
        elif(xx.name == 's1'):
            LS1 += [xx]
        elif(xx.name == 's2'):
            LS2 += [xx]
        elif(xx.name == 'w'):
            LW += [xx]
        elif(xx.name == 't'):          # This will not appear in the actual process, we add it here to facilitate debugging
            LT += [xx]
        else:
            print "Problem in permu_1_va!!!\n\n\n"
            pdb.set_trace()
    AA = permu_general(LT1)
    BB = permu_general(LT2)
    BB3 = permu_general(LT3)
    CC = permu_general(LS1)
    DD = permu_general(LS2)
    EE = permu_general(LW)
    FF = permu_general(LDEL)
    GG = permu_general(LT)
    ZZ = []
    for x in [AA, BB, BB3, CC, DD, EE, FF, GG]:
        if(len(x) != 0):
            ZZ += [x]
    LLL = combi_lists(ZZ) # LLL consists of lists
    for xxx in LLL:
        yyy = []
        for xyx in xxx:
            yyy += xyx
        mats = LV + LDEN + yyy
        result += [Term(Coefficient(numm, mats), ops)]
    return result



def unique_section_divide(cons, section_size, section_num):

    """ return a list, whose element is a section_num of sections which contain section_size elements from cons """
    """ for element in one section, permutation is not considered """
    """ for every element in the final list, permutation between sections are not considered """
    """ so we need to find a unique way to divide """

    #pdb.set_trace()
    result = []
    if(len(cons) < section_size):
        result = []
    elif(section_num == 0):                                                                        #  stop RECURSION
        result = []
    elif(section_num == 1):
        result_prep =  combi(section_size, cons)
        for xx in result_prep:
            result += [ [xx] ]                                                                     # format at any recursive step must be the same as the final result
    else:
        L1 = combi(section_size, cons)                                                             # start RECURSION
        for xx in L1:
            cons_left = []
            for yy in cons:
                if(yy not in xx and yy > xx[0]):                                                   # we require the first element in next section increase to achieve
                    cons_left += [yy]                                                              # uniqueness
            the_other_parts = unique_section_divide(cons_left, section_size, section_num-1)        # RECURSION point
            for zz in the_other_parts:
                result += [ [xx] + zz ]
    return result


def permu_unique_section_divide(cons, section_size, section_num):

    """ cons = [2, 4, 5, 6, 7], section_size = 2, section_num = 2, return [      [ [2,5], [6,7] ], ...   ]"""
    """ return a list, whose element is a list, which contains section_num of sections which contain section_size elements from cons """
    """ for every element in the final list, permutation between sections are allowd """
    """ for element in one section, permutation is not considered """

    debug = 0
    if(debug):
        print  "\n debug permu_unique_section_divide \n"
        pdb.set_trace()
    result = []
    if(section_num == 0):                                                                          #  stop RECURSION
        result = []
    elif(section_num == 1):
        result_prep =  combi(section_size, cons)
        for xx in result_prep:
            result += [ [xx] ]
    else:
        L1 = combi(section_size, cons)                                                             # start RECURSION
        for xx in L1:
            cons_left = deletes(cons, give_offset_s(xx, cons))
            the_other_parts = permu_unique_section_divide(cons_left, section_size, section_num-1)  # RECURSION point
            for zz in the_other_parts:
                result += [ [xx] + zz ]
    return result


def inter_out_permu_unique_section_divide(cons, section_size, section_num):                         # outside permu followed by internal permu

    """ cons = [2, 4, 5, 6, 7], section_size = 2, section_num = 2, return [      [ [2,5], [6,7] ], ...   ]"""
    """ return a list, whose element is a list, which contains section_num of sections which contain section_size elements from cons """
    """ for every element in the final list, permutation between sections are allowd """
    """ for element in one section, internal permutation is considered """

    debug = 0
    if(debug):
        print  "\n debug permu_unique_section_divide \n"
        pdb.set_trace()
    result = []
    if(len(cons) < section_size):
        result = []
    elif(section_num == 0):                                                                          #  stop RECURSION
        result = []
    elif(section_num == 1):
        result_prep =  permu(section_size, cons)
        for xx in result_prep:
            result += [ [xx] ]
    else:
        L1 = permu(section_size, cons)                                                               # start RECURSION
        for xx in L1:
            cons_left = deletes(cons, give_offset_s(xx, cons))
            the_other_parts = inter_out_permu_unique_section_divide(cons_left, section_size, section_num-1)  # RECURSION point
            for zz in the_other_parts:
                result += [ [xx] + zz ]
    return result



def inter_out_permu_unique_section_divide_old(cons, section_size, section_num):                      # outside permu followed by internal permu

    """ cons = [2, 4, 5, 6, 7], section_size = 2, section_num = 2, return [      [ [2,5], [6,7] ], ...   ]"""
    """ return a list, whose element is a list, which contains section_num of sections which contain section_size elements from cons """
    """ for every element in the final list, permutation between sections are allowd """
    """ for element in one section, internal permutation is considered """

    step1 = permu_unique_section_divide(cons, section_size, section_num)
    result = []
    for xx in step1:
        result += inter_permu_list(xx)
    return result



def inter_permu_list(one_con):

    """ one_con = [ [2,5], [6,7] ], return [     [ [2,5], [6,7]], [[5,2], [6,7]],  [[2,5], [7,6]],  [[5,2], [7,6]]    ]  """

    size = len(one_con[0])
    total = []
    for xx in one_con:
        total += [permu(size, xx)]
    return combi_pick_one(total)


def combi_pick_one(listss):

    """ pick one from every list in lists to form a list; return all lists from this process """
    """ listss = [ [[2,5], [5, 2]], [[6,7], [7,6]] ], return  [ [[2,5], [6,7]], [...], [...], [...]  ]"""

    result = []
    leng = len(listss)
    if(leng == 1):
        for xx in listss[0]:
            result += [ [xx] ]
    elif(leng == 0):
        pdb.set_trace()
    else:
        for xx in listss[0]:
            for yy in combi_pick_one(listss[1:]):
                result += [ [xx] + yy ]
    return result



    


    

    
  






#**************************************************
#               symmetry_number                   *
#**************************************************





def symmetry_number(gterm):                                     # gives the symmetry number n! (float) for gterm depending on the permutation
    N = len(gterm.uOperators[0].upperIndicees)
    if(N in [0, 1]):
        return 1.0
    elif(N in [2, 3, 4, 5]):
        return symmetry_number_help(gterm)
    else:
        print "\n\n                      Something unassigned in symmetry_number@real.py \n\n"
        printterm(gterm)
        pdb.set_trace()
        X = im_triple_exclude_G_1(gterm)
        print X



def symmetry_number_help(gterm):
    N = len(gterm.uOperators[0].upperIndicees)
    if(N == 2):
        coll = [[gterm.uOperators[0].upperIndicees[0].type, gterm.uOperators[0].upperIndicees[0].active,
                 gterm.uOperators[0].lowerIndicees[0].type, gterm.uOperators[0].lowerIndicees[0].active],
                [gterm.uOperators[0].upperIndicees[1].type, gterm.uOperators[0].upperIndicees[1].active,
                 gterm.uOperators[0].lowerIndicees[1].type, gterm.uOperators[0].lowerIndicees[1].active]]
    elif(N == 3):
        coll = [[gterm.uOperators[0].upperIndicees[0].type, gterm.uOperators[0].upperIndicees[0].active,
                 gterm.uOperators[0].lowerIndicees[0].type, gterm.uOperators[0].lowerIndicees[0].active],
                [gterm.uOperators[0].upperIndicees[1].type, gterm.uOperators[0].upperIndicees[1].active,
                 gterm.uOperators[0].lowerIndicees[1].type, gterm.uOperators[0].lowerIndicees[1].active],
                [gterm.uOperators[0].upperIndicees[2].type, gterm.uOperators[0].upperIndicees[2].active,
                 gterm.uOperators[0].lowerIndicees[2].type, gterm.uOperators[0].lowerIndicees[2].active]]
    elif(N == 4):
        coll = [[gterm.uOperators[0].upperIndicees[0].type, gterm.uOperators[0].upperIndicees[0].active,
                 gterm.uOperators[0].lowerIndicees[0].type, gterm.uOperators[0].lowerIndicees[0].active],
                [gterm.uOperators[0].upperIndicees[1].type, gterm.uOperators[0].upperIndicees[1].active,
                 gterm.uOperators[0].lowerIndicees[1].type, gterm.uOperators[0].lowerIndicees[1].active],
                [gterm.uOperators[0].upperIndicees[2].type, gterm.uOperators[0].upperIndicees[2].active,
                 gterm.uOperators[0].lowerIndicees[2].type, gterm.uOperators[0].lowerIndicees[2].active],
                [gterm.uOperators[0].upperIndicees[3].type, gterm.uOperators[0].upperIndicees[3].active,
                 gterm.uOperators[0].lowerIndicees[3].type, gterm.uOperators[0].lowerIndicees[3].active]]
    elif(N == 5):
        coll = [[gterm.uOperators[0].upperIndicees[0].type, gterm.uOperators[0].upperIndicees[0].active,
                 gterm.uOperators[0].lowerIndicees[0].type, gterm.uOperators[0].lowerIndicees[0].active],
                [gterm.uOperators[0].upperIndicees[1].type, gterm.uOperators[0].upperIndicees[1].active,
                 gterm.uOperators[0].lowerIndicees[1].type, gterm.uOperators[0].lowerIndicees[1].active],
                [gterm.uOperators[0].upperIndicees[2].type, gterm.uOperators[0].upperIndicees[2].active,
                 gterm.uOperators[0].lowerIndicees[2].type, gterm.uOperators[0].lowerIndicees[2].active],
                [gterm.uOperators[0].upperIndicees[3].type, gterm.uOperators[0].upperIndicees[3].active,
                 gterm.uOperators[0].lowerIndicees[3].type, gterm.uOperators[0].lowerIndicees[3].active],
                [gterm.uOperators[0].upperIndicees[4].type, gterm.uOperators[0].upperIndicees[4].active,
                 gterm.uOperators[0].lowerIndicees[4].type, gterm.uOperators[0].lowerIndicees[4].active]]
    else:
        printterm(gterm)
        pdb.set_trace()
    permucoll = permu(N, coll)
    symmetrynumber = 0
    for xx in permucoll:
        if(xx == coll):
            symmetrynumber += 1
    return float(symmetrynumber)



def coeff_from_mats(term, names):                # determine the coefficient.const from the component mats.
    result = 1
    matname = []
    mats_time = []
    for xx in term.coefficient.matElement:
        if(xx.name in matname):
            mats_time[matname.index(xx.name)] += 1
        elif(xx.name in names):
            matname += [xx.name]
            mats_time += [1]
    for yy in mats_time:
        result *= 1.0/factorial(yy)
    return result





#**************************************************
#                   specials                      *
#**************************************************





def setchar(term):                     # set attribute to a term
    if(len(term.uOperators) == 1):
        if(len(term.uOperators[0].upperIndicees) == 1):
            kk = term.uOperators[0].upperIndicees[0]
            gg = term.uOperators[0].lowerIndicees[0]
            if(kk.type == "p" and kk.num == "1" and gg.type == "h" and gg.num == "1"):
                term.uOperators[0].cha = 1
            elif(kk.type == "p" and kk.num == "1" and gg.type == "e" and gg.num == "1"):
                term.uOperators[0].cha = 2
            elif(kk.type == "m" and kk.num == "1" and gg.type == "h" and gg.num == "1"):
                term.uOperators[0].cha = 3
        elif(len(term.uOperators[0].upperIndicees) == 0):
            term.uOperators[0].cha = 0
        elif(len(term.uOperators[0].upperIndicees) == 2):
            kk1 = term.uOperators[0].upperIndicees[0]
            kk2 = term.uOperators[0].upperIndicees[1]
            gg1 = term.uOperators[0].lowerIndicees[0]
            gg2 = term.uOperators[0].lowerIndicees[1]
            if(kk1.type == "p" and kk1.num == "1" and gg1.type == "h" and gg1.num == "1" and kk2.type == "p" and kk2.num == "2" and gg2.type == "h" and gg2.num == "2"):
                term.uOperators[0].cha = 4
            elif(kk1.type == "p" and kk1.num == "1" and gg1.type == "e" and gg1.num == "1" and kk2.type == "p" and kk2.num == "2" and gg2.type == "h" and gg2.num == "1"):
                term.uOperators[0].cha = 5
            elif(kk1.type == "p" and kk1.num == "1" and gg1.type == "h" and gg1.num == "1" and kk2.type == "m" and kk2.num == "1" and gg2.type == "h" and gg2.num == "2"):
                term.uOperators[0].cha = 6
    return term





def setchars(terms):                   # set attribute to each term
    res = []
    for xx in terms:
        res = res + [setchar(xx)]
    return res



def pickupchar(terms, ch):             # pick up some terms of particular attribute
    tt = setchars(terms)
    res = []
    for xx in tt:
        if(len(xx.uOperators) == 1):
            if(xx.uOperators[0].cha == ch):
                res = res + [xx]
    return res



def set(term1, upp, low):              

    """ for terms containing generic indices, set them to be corresponding terms specified by upp and low, upp/low is, e.g., [1, 0, 0, 1, 0] """
    """ I assume this function is not in use and forgot how it works and what it does """

    term = copy.deepcopy(term1)       
    if(len(choosetype(term.uOperators[0].upperIndicees, "g") + choosetype(term.uOperators[0].lowerIndicees, "g")) == 0):
        return term # if there is no 'g' indices, return term
    else:
        list = []
        cxa = counttype(term.uOperators[0].upperIndicees)
        cxb = upp
        for aa in range(4):
              if(cxa[aa] > cxb[aa]):
                  return term          # compare upperindicees with upp, if term.uOperators[0].upperindices contains more indices of any kind than upp, 
                                       # return, because this term then cann't belong to the class speicified by upp and low
        cya = counttype(term.uOperators[0].lowerIndicees)
        cyb = low
        for aa in range(4):
              if(cya[aa] > cyb[aa]):
                  return term          # do the same to lowerindices as above
        typelist = ["p", "e", "h", "m", "g"]
        nt = numtype(term)
        cx = []
        for bb in range(4):
            cx = cx + [cxb[bb]-cxa[bb]]
        for cc in range(4):
            for dd in range(cx[cc]):
                list = list + [Index(typelist[cc], str(nt[cc]+1+dd))]
        cy = []
        for ee in range(4):
            cy = cy + [cyb[ee]-cya[ee]]
        for ff in range(4):
            for gg in range(cy[ff]):
                list = list + [Index(typelist[ff], str(nt[ff]+gg+1))]
        glist = choosetype(term.uOperators[0].upperIndicees, "g") + choosetype(term.uOperators[0].lowerIndicees, "g")
        ind = makeIndTerm(term)
        for a in range(len(ind.indexSets)):
            for b in range(len(ind.indexSets[a].indexPairs)):
                if(ind.indexSets[a].indexPairs[b].upperIndex.type == "g"):
                    ind.indexSets[a].indexPairs[b].upperIndex = list[ind.indexSets[a].indexPairs[b].upperIndex.gvIndexNo(glist)]
                if(ind.indexSets[a].indexPairs[b].lowerIndex.type == "g"):
                    ind.indexSets[a].indexPairs[b].lowerIndex = list[ind.indexSets[a].indexPairs[b].lowerIndex.gvIndexNo(glist)]
        term = arrangeIndicees(ind, term)
        return term




def setterms(terms, upp, low):         # set function acting on terms
    result = []
    for xx in terms:
        result = result + [set(copy.deepcopy(xx), upp, low)]
    return result



def setfree(term1, list):              # set indices in 'list' free
    term = copy.deepcopy(term1)
    ind = makeIndTerm(term)
    for aa in ind.indexSets:
        for bb in aa.indexPairs:
            if(bb.upperIndex.inArray(list)):
                bb.upperIndex.fix = 0
            if(bb.lowerIndex.inArray(list)):
                bb.lowerIndex.fix = 0
    term = arrangeIndicees(ind, term)
    return term



def setfrees(terms, list):            
    res = []
    for xx in terms:
        res = res + [setfree(xx,list)]
    return res


def setfix(term1, list):               # set some indices to be fixed
    term = copy.deepcopy(term1)
    ind = makeIndTerm(term)
    for aa in ind.indexSets:
        for bb in aa.indexPairs:
            if(bb.upperIndex.inArray(list)):
                bb.upperIndex.fix = 1
            if(bb.lowerIndex.inArray(list)):
                bb.lowerIndex.fix = 1
    term = arrangeIndicees(ind, term)
    return term




def setfixs(terms, list):            
    res = []
    for xx in terms:
        res = res + [setfix(xx,list)]
    return res



def reset(term1, pre, aft):            # Suppose there are overlap between pre and aft, say pre[1] == aft[0]; if there are repeated indices in term1, say a = b = pre[0], after a is converted
    """ replace "pre" indices with "aft" indices """  # to aft[0], i guess b is also converted to aft[0] automatically (?) if they are references to each other. then when it comes to 'b',
    if(len(pre) == 0):                                # it has value aft[0]; since now aft[0] belongs to 'pre', shouldn't it be converted again to aft[1]? If this is true, this is wrong.
        return term1                                  # but so far no problem aries yet....
    else:                                             # Now I understand: reassignment will cause the interpreter to allocate a new address for the lv, but the other varialbes pointing to
        term = copy.deepcopy(term1)                   # same address are not affected. Reassignment is not an operation on the address. 
        ind = makeIndTerm(term)
        for a in range(len(ind.indexSets)):
            for b in range(len(ind.indexSets[a].indexPairs)):
                    if(ind.indexSets[a].indexPairs[b].upperIndex.inArray(pre)):
                        ind.indexSets[a].indexPairs[b].upperIndex = aft[ind.indexSets[a].indexPairs[b].upperIndex.gvIndexNo(pre)]
                    if(ind.indexSets[a].indexPairs[b].lowerIndex.inArray(pre)):
                        ind.indexSets[a].indexPairs[b].lowerIndex = aft[ind.indexSets[a].indexPairs[b].lowerIndex.gvIndexNo(pre)]
        term = arrangeIndicees(ind, term)
        return term






def resets(terms, pre, aft):           # replace "pre" indices with "aft" indices
    res = []
    for xx in terms:
        res = res + [reset(xx, pre, aft)]
    return res



def inc_index(term1, type, step):      # increase index.num by step for type "type", to avoid overlapping naming
    term = copy.deepcopy(term1)
    us = term.uOperators[0].upperIndicees
    ls = term.uOperators[0].lowerIndicees
    for xx in range(len(us)):
        if us[xx].fix != 1:
            if us[xx].type == type: us[xx] = Index(type, str(int(us[xx].num) + step), us[xx].att, us[xx].fix, us[xx].active, us[xx].other)
        if ls[xx].fix != 1:
            if ls[xx].type == type: ls[xx] = Index(type, str(int(ls[xx].num) + step), ls[xx].att, ls[xx].fix, ls[xx].active, ls[xx].other)
    for xx in term.coefficient.matElement:
        xxus = xx.matUpperIndicees
        xxls = xx.matLowerIndicees
        for yy in range(len(xx.matUpperIndicees)):
            if xxus[yy].fix != 1:
                if xxus[yy].type == type: xxus[yy] = Index(type, str(int(xxus[yy].num) + step), xxus[yy].att, xxus[yy].fix, xxus[yy].active, xxus[yy].other)
            if xxls[yy].fix != 1:
                if xxls[yy].type == type: xxls[yy] = Index(type, str(int(xxls[yy].num) + step), xxls[yy].att, xxls[yy].fix, xxls[yy].active, xxls[yy].other)
    return term





def inc_index_s(terms, type, step):    # increase index.num by step for type "type", to avoid overlapping naming
    res = []
    for xx in terms: res.append(inc_index(xx, type, step))
    return res



# expand specific particle indices into 'g' and 'h'
def expand_specific_p_1(term1, p1):      
    term = copy.deepcopy(term1)   
    indlist = give_indlist(term)
    h = 0
    g = 0
    for xx in indlist:
        if(xx.type == 'h' and int(xx.num) > h):
            h = int(xx.num)
        elif(xx.type == 'g' and int(xx.num) > g):
            g = int(xx.num)
    newg = Index("g", str(g+1))
    newh = Index("h", str(h+1))
    return [reset(term, [p1], [newg])] + multiply_co(-1.0, [reset(term, [p1], [newh])])


# expand specific particle indices into 'g' and 'h'
def expand_s_specific_p_1(term1s, p1):      
    res = []
    for xx in term1s:
        res += expand_specific_p_1(xx, p1)
    return res



def expand_specific_ps(term1, ps):      
    term = copy.deepcopy(term1)   
    res = [term]
    for xx in ps:
        res = expand_s_specific_p_1(res, xx)
    return res




def expand_specific_g(term1, gs):      
    term = copy.deepcopy(term1)   
    indlist = give_indlist(term)
    p = 0
    h = 0
    g = 0
    for xx in indlist:
        if(xx.type == 'p' and int(xx.num) > p):
            p = int(xx.num)
        elif(xx.type == 'h' and int(xx.num) > h):
            h = int(xx.num)
    if len(gs) == 0:
        return [term]
    else:
        glist = gs
        pl = []
        hl = []
        for i in range(len(gs)):
            pl += [Index("p", str(p+1+i))]
            hl += [Index("h", str(h+1+i))]
        if(len(glist) == 1):
            r_list = [[pl[0]],
                      [hl[0]]]
        elif(len(glist) == 2):
            r_list = [ [pl[0], pl[1]],
                       [hl[0], hl[1]],
                       [pl[0], hl[0]],
                       [hl[0], pl[0]] ]
        elif(len(glist) == 3):
            r_list = [ pl[:3],
                       [pl[0], pl[1], hl[0]],
                       [pl[0], hl[0], pl[1]],
                       [hl[0], pl[0], pl[1]],
                       [pl[0], hl[0], hl[1]],
                       [hl[0], pl[0], hl[1]],
                       [hl[0], hl[1], pl[0]],
                       hl[:3] ]
        elif(len(glist) == 4):
            r_list = [ pl,
                       pl[:3] + [hl[0]],
                       [pl[0], pl[1], hl[0], pl[2]],
                       [pl[0], hl[0], pl[1], pl[2]],
                       [hl[0]] + pl[:3],
                       [pl[0], pl[1], hl[0], hl[1]],
                       [pl[0], hl[0], pl[1], hl[1]],
                       [pl[0], hl[0], hl[1], pl[2]],
                       [hl[0], pl[0], pl[1], hl[1]],
                       [hl[0], pl[0], hl[1], pl[1]],
                       [hl[0], hl[1], pl[0], pl[1]],
                       hl[:3] + [pl[0]],
                       hl[:2] + [pl[0], hl[2]],
                       [hl[0], pl[0], hl[1], hl[2]],
                       [pl[0]] + hl[:3],
                       hl ]
        else:
            print len(glist)
            print "\n error in break_g@eq_lib\n"
            pdb.set_trace()
        res = []
        for xx in r_list:
            res += [reset(term, glist, xx)]
        return res





def expand_specific_g_s(terms, gs):             
    res = []
    for xx in terms:
        res += expand_specific_g(xx, gs)
    return res




def break_g(term1):                    # convert "g" indices to all possible p/h('e' or 'm' is included in p/h) indices,  
    term = copy.deepcopy(term1)        # this method may be used to replace set(), setterms() fucntions
    indlist = give_indlist(term)
    p = 0
    h = 0
    g = 0
    for xx in indlist:
        if(xx.type == 'p' and int(xx.num) > p):
            p = int(xx.num)
        elif(xx.type == 'h' and int(xx.num) > h):
            h = int(xx.num)
        elif(xx.type == 'g'):
            g += 1
    if(g == 0):
        return [term]
    else:
        #glist = choosetype(term.uOperators[0].upperIndicees, "g") + choosetype(term.uOperators[0].lowerIndicees, "g")
        glist = choosetype(indlist, "g") #  the above line is for traditional cases; once going to km_gwt, that is not enough
        pl = []
        hl = []
        for i in range(4):
            pl += [Index("p", str(p+1+i))]
            hl += [Index("h", str(h+1+i))]
        if(len(glist) == 1):
            r_list = [[pl[0]],
                      [hl[0]]]
        elif(len(glist) == 2):
            r_list = [ [pl[0], pl[1]],
                       [hl[0], hl[1]],
                       [pl[0], hl[0]],
                       [hl[0], pl[0]] ]
        elif(len(glist) == 3):
            r_list = [ pl[:3],
                       [pl[0], pl[1], hl[0]],
                       [pl[0], hl[0], pl[1]],
                       [hl[0], pl[0], pl[1]],
                       [pl[0], hl[0], hl[1]],
                       [hl[0], pl[0], hl[1]],
                       [hl[0], hl[1], pl[0]],
                       hl[:3] ]
        elif(len(glist) == 4):
            r_list = [ pl,
                       pl[:3] + [hl[0]],
                       [pl[0], pl[1], hl[0], pl[2]],
                       [pl[0], hl[0], pl[1], pl[2]],
                       [hl[0]] + pl[:3],
                       [pl[0], pl[1], hl[0], hl[1]],
                       [pl[0], hl[0], pl[1], hl[1]],
                       [pl[0], hl[0], hl[1], pl[2]],
                       [hl[0], pl[0], pl[1], hl[1]],
                       [hl[0], pl[0], hl[1], pl[1]],
                       [hl[0], hl[1], pl[0], pl[1]],
                       hl[:3] + [pl[0]],
                       hl[:2] + [pl[0], hl[2]],
                       [hl[0], pl[0], hl[1], hl[2]],
                       [pl[0]] + hl[:3],
                       hl ]
        else:
            print len(glist)
            print "\n error in break_g@eq_lib\n"
            pdb.set_trace()
        res = []
        for xx in r_list:
            res += [reset(term, glist, xx)]
        return res





# originally named break_g_s()
def expand_general_indices(terms):             
    res = []
    for xx in terms:
        res += break_g(xx)
    return res



def set_inactive(term1, index):        # set some index in term to be inactive, also refer to texform() to see the difference of output
    term = copy.deepcopy(term1)
    ind = makeIndTerm(term)
    for xx in ind.indexSets:
        for yy in xx.indexPairs:
            if(yy.upperIndex.equal(index)):
                yy.upperIndex.active = -1
            if(yy.lowerIndex.equal(index)):
                yy.lowerIndex.active = -1
    term = arrangeIndicees(ind, term)
    return term




def set_inactive_s(terms, index):
    res = []
    for x in terms:
        res += [set_inactive(x, index)]
    return res




def set_term(term1):
    term = copy.deepcopy(term1)
    type = []
    if(len(term.uOperators[0].upperIndicees) == 0):
        term.type = 0
        return term
    else:
        for aa in term.uOperators[0].upperIndicees:
            type += [aa.type]
        for bb in term.uOperators[0].lowerIndicees:
            type += [bb.type]
        for ab in a_type:
            if(type == ab):
                term.type = a_type.index(ab)
                break
        return term





def set_terms(terms):                  # set the value of term.type for terms
    res = []
    for xx in terms:
        res += [set_term(xx)]
    return res


def D_2_OP(terms):

    """ convert the density matrix to the operator """

    res = []
    for xx in terms:
        res += [D_2_OP_1(xx)]
    return res


def D_2_OP_1(term1):

    """ convert the density matrix to the operator """

    term = copy.deepcopy(term1)
    A = []
    B = []
    for xx in term.coefficient.matElement:
        if(xx.name in ['d1', 'd2', 'd3', 'd4']):
            A += [xx]
        else:
            B += [xx]
    if(len(A) == 0):
        return term
    else:
        return Term(Coefficient(term.coefficient.const, B), [Uoperator(A[0].matUpperIndicees, A[0].matLowerIndicees)], term.type, term.other)



def delete_x_s(terms):
    R = []
    for xx in terms:
        R += [delete_x(xx)]
    return R


def delete_x(term):                    # 'X' is a hypothetical mat which helps to reduce # of terms in Derivative temparily plugged in, 
    R = []                             # so we must delete it finally
    for xx in term.coefficient.matElement:
        if(xx.name != 'X'):
            R += [xx]
    return Term(Coefficient(term.coefficient.const, R), term.uOperators, term.type, term.other)  





#**************************************************
#             product and commutator             *
#**************************************************





def formproduct(term1, term2):

        """ return the product of term1 and term2 """

        if 0:
            set_trace()
            printterms([term1, term2])
        term = Term(Coefficient(1.0, []), [Uoperator([], [])])
        terma = copy.deepcopy(term1)
        termb = copy.deepcopy(term2)
        termb = inc_index_avoid_overlap([termb], [terma])[0]
        terma = inc_index_avoid_overlap([terma], [termb])[0]                   # sometimes 'termb' contains fixed indices; mostly 'terma' does
        term.coefficient.const = terma.coefficient.const * termb.coefficient.const
        term.coefficient.matElement = terma.coefficient.matElement + termb.coefficient.matElement
        if(len(term1.uOperators[0].upperIndicees) == 0 and len(term2.uOperators[0].upperIndicees) == 0): term.uOperators = [Uoperator([], [])]
        elif(len(term1.uOperators[0].upperIndicees) == 0): term.uOperators = termb.uOperators
        elif(len(term2.uOperators[0].upperIndicees) == 0): term.uOperators = terma.uOperators
        else:
                term.uOperators = terma.uOperators + termb.uOperators
        return term





def product(A, B):                     # A, B are terms, return the product of terms in A and B
    res = []
    for a in A:
            print "a = "
            printterms([a])
            for b in B:
                    #print "b = "
                    res += [copy.deepcopy(formproduct(a, b))]
                    #print "result = "
                    #printterms(res)
    return res


def normal_product_s(termslist2):
    res = []
    for a in termslist2[0]:
        for b in termslist2[1]:
            res += [normal_2product([a, b])]
    return res
    


def normal_2product(terrms):

    """ return {AB}, {} denotes normal order """
    terma = copy.deepcopy(terrms[0])
    termb = copy.deepcopy(terrms[1])
    termb = inc_index_avoid_overlap([termb], [terma])[0]
    terma.coefficient.const *= termb.coefficient.const
    terma.coefficient.matElement += termb.coefficient.matElement
    terma.uOperators[0].upperIndicees += termb.uOperators[0].upperIndicees
    terma.uOperators[0].lowerIndicees += termb.uOperators[0].lowerIndicees
    return terma





def normal_product(terrms):

    """ return {ABC...}, {} denotes normal order """

    term = Term(Coefficient(1.0, []), [Uoperator([], [])])
    terms = copy.deepcopy(terrms)
    term.coefficient.const = 1.0
    for xx in terms:
        term.coefficient.const *= xx.coefficient.const
        term.coefficient.matElement += xx.coefficient.matElement
        term.uOperators[0].upperIndicees += xx.uOperators[0].upperIndicees
        term.uOperators[0].lowerIndicees += xx.uOperators[0].lowerIndicees
    return term



def normal_order_term_term(ta, tb):       # write ta \times tb in normal order
    if(len(ta.uOperators[0].upperIndicees) == 0):
        return [Term(Coefficient(ta.coefficient.const * tb.coefficient.const, ta.coefficient.matElement + tb.coefficient.matElement), tb.uOperators)] 
    elif(len(tb.uOperators[0].upperIndicees) == 0):
        return [Term(Coefficient(ta.coefficient.const * tb.coefficient.const, ta.coefficient.matElement + tb.coefficient.matElement), ta.uOperators)]
    else:
        gwc = formproduct(ta, tb)
        theochem1 = solvetwoop(gwc)
        return theochem1



def normal_order_term_terms(ta, Ts):    
    result = []
    for xx in Ts:
        result += normal_order_term_term(ta,xx)
    return result



def normal_order_terms_term(Ts, ta):    
    result = []
    for xx in Ts:
        result += normal_order_term_term(xx, ta)
    return result



def normal_order_terms_terms(Ts1, Ts2): 
    result = []
    for xx in Ts1:
        result += normal_order_term_terms(xx, Ts2)
    return result



def contraction_term_term(ta, tb):       

    """ evaluate (a b)_connected, a and b are two terms, return a list of terms """

    if(len(ta.uOperators[0].upperIndicees) == 0 or len(tb.uOperators[0].upperIndicees) == 0):
        print "                        one of the two terms is a constant in contraction_term_term@real.py"
        printterms([ta, tb])
        return []
    else:
        gwc = formproduct(ta, tb)
        lengg = len(ta.uOperators[0].upperIndicees) + len(tb.uOperators[0].upperIndicees)
        theochem1 = solvetwoop(gwc)
        theochem2 = []
        for jcp in theochem1:
            if(len(jcp.uOperators[0].upperIndicees) < lengg):
                theochem2 += [jcp]
        #return merge_all(theochem2)   # G have 3-body components, merge_all may not work
        return theochem2



def contraction_term_terms(ta, Ts):    # evaluate (a A)_c
    result = []
    clock = 1
    for xx in Ts:
        #print clock
        result += contraction_term_term(ta,xx)
        clock += 1
    return result



def contraction_terms_term(Ts, ta):    # evaluate (A a)_c
    result = []
    for xx in Ts:
        result += contraction_term_term(xx, ta)
    return result



def contraction_terms_terms(Ts1, Ts2): # evaluate (A B)_c
    result = []
    #print "total # of terms: to be processed@contraction_terms_terms: " + str(len(Ts1) * len(Ts2)) + '\n\n'
    #counter = 0
    for xx in Ts1:
        result += contraction_term_terms(xx, Ts2)
        #if counter % 10 == 0 or 1: print "the " + str(counter) + "-th term\n"
        #counter += 1
    #print "done\n\n"
    return result





def commutator(A, B):
    """  return [A, B] """

    print "\n                    @@@@@@@@                commutator starts                @@@@@@@@\n"
    res = []
    if len(A) == 0 or len(B) == 0: return res
    else:
        first = []
        second = []
        a = removezero(A)
        b = removezero(B)
        i = 0                                                                
        for ab in a:                                                    
            for bc in b:
                i += 1
                if i - 100 * (i/100) == 0: print "       J  =  %s    "   %   i
                ra = len(ab.uOperators[0].upperIndicees)
                rb = len(bc.uOperators[0].upperIndicees)
                first += filter(solvetwoop(formproduct(ab, bc)), ra + rb - 1)
                second += filter(solvetwoop(formproduct(bc, ab)), ra + rb - 1)
        for xx in first: res.append(copy.deepcopy(xx))
        if(len(second) == 0): pass
        else:
            for i in range(len(second)):
                cd = copy.deepcopy(second[i])
                cd.coefficient.const *= -1.0
                res.append(cd)
        print "\n                    @@@@@@@@                commutator ends                  @@@@@@@@\n"
        return res









#**************************************************
#                  selecting                      *
#**************************************************



def discard_3d_term_s(terms):                                                 # discard terms containing 3 particle DM
    res = []
    for xx in terms:
        if(contain_3dm(xx) == 1):
            pass
        else:
            res += [xx]
    return res


def contain_3dm(term):                                                        # determine if there is a 3DM in mats
    res = 0
    for xx in term.coefficient.matElement:
        if(xx.name in ['nd', 'd1', 'd2', 'd3', 'd4']):
            if(len(xx.matUpperIndicees) == 3):
                res = 1
                break
    return res




def discard_terms_have_mat(terms, matname):                                    # discard terms which have some matelement
    res = []
    for xx in terms:
        flag = 1
        for yy in xx.coefficient.matElement:
            if(yy.name == matname):
                flag = 0
                break
        if(flag == 1):
            res += [xx]
    return res
   

def select_terms_rank_n(terms, namerange1, namerange2,  n=0):                    
    """ select terms of rank n with respect to matrixelements of names in 'namerange1', say. ['s1', 's3'], while not containing terms of name 
        in 'namerange2' """
    
    res = []                                                        
    for term in terms:                                             
        num1 = 0                                                  
        num2 = 0
        for xx in term.coefficient.matElement:
            if xx.name in namerange1: num1 += 1
            elif xx.name in namerange2: num2 += 1
        if num1 == n and num2 == 0: res.append(term)
        elif len(namerange1) == 0 and num2 == 0:                              # this is to discard the terms with mat name in 'namerange2'
            res.append(term)
    return res  


def pick_op_rank(terms, N):
    return [xx for xx in terms if len(xx.uOperators[0].upperIndicees) == N]




def choosetype(indlist, x):                                                    # return indices in list of type x
    res = []
    for i in range(len(indlist)):
        if(indlist[i].type == x):
            res = res + [indlist[i]]
    return res




def give_type(term):                                                           # give the list of type of indices of operator of term
    result = ""
    for aa in term.uOperators[0].upperIndicees:
        result += aa.type
    for bb in term.uOperators[0].lowerIndicees:
        result += bb.type
    return result



def give_mattype(term):
    result = []
    for aa in term.coefficient.matElement[0].matLowerIndicees:
        result += [aa.type]
    for bb in term.coefficient.matElement[0].matUpperIndicees:
        result += [bb.type]
    return result



def pick_type(terms, type):                                                    # pick up terms of type "type"
    result = []
    for xx in terms:
        if(give_type(xx) == type):
            result += [xx]
    return result



def choose_mat_name(term, name1):                                              # obtain matrixelemens of the name "name1"
    result = []
    for x in term.coefficient.matElement:
        if(x.name == name1):
             result += [copy.deepcopy(x)]
    return result    

def choose_terms_matname(terms, namee):                                        # select terms containing an matrixelement named 'name'
    res = []
    for xx in terms:
        flag = 0
        for yy in xx.coefficient.matElement:
            if(yy.name == namee):
                flag = 1
                break
        if(flag):
            res += [xx]
    return res



def choose_terms_matnames(terms, namees):                                        # select terms containing an matrixelement named 'name'
    res = []
    for xx in terms:
        flag = 0
        for yy in xx.coefficient.matElement:
            if(yy.name in namees):
                flag = 1
                break
        if flag: res.append(xx)
    return res





def filter_type(terms, types, n):       # discard terms whose operators have more than n 'type' of indices
    result = []
    for xx in terms:
        if filter_type_judge(xx, types, n) == 0: result.append(xx)
    return result


def filter_type_judge(term, types, n):
    nz = counttype(term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees)
    if(len(types) == 1):
        if(nz[['p','e', 'h', 'm', 'g'].index(type)] > n):
            return 1
        else:
            return 0
    elif(len(types) == 2):
        if(nz[['p','e', 'h', 'm', 'g'].index(types[0])] + nz[['p','e', 'h', 'm', 'g'].index(types[1])] > n):
            return 1
        else:
            return 0
    else:
        pdb.set_trace()
    
     



def independent_mat_help(mats1, mats2):  # divide mats2 to 2 classes according to whether being contracted with mats1 or not, return contracted ones and
    NY = []                              # uncontacted ones
    WS = []
    if(len(mats2) == 0):
        return [[], []]
    else:
        ind = []
        for x in mats1:
            ind += x.matUpperIndicees + x.matLowerIndicees
        for xx in mats2:
            flag = 0
            ind2 = xx.matUpperIndicees + xx.matLowerIndicees
            for yy in range(len(ind2)):
                if(ind2[yy].inArray(ind)):
                    flag = 1
                    break
            if(flag):
                NY += [xx]
            else:
                WS += [xx]
        return [NY, WS]



def independent_mat(term1):             # find the uncontracted matrix element
    term = copy.deepcopy(term1)
    PR = []
    UCB = []
    find_h = 0
    for xx in term.coefficient.matElement:
        if(xx.name in ['f', 'v', 'nf', 'nv', 'kmf', 'kmv', 'hbar1', 'hbar2', 'lvsA', 'lvsB', 'lvsC', 'lvsD'] and system_switch['ccsd(t)'] == 0):
            PR += [xx]
        elif(xx.name in ['f', 'v', 'nf', 'nv', 'kmf', 'kmv', 'hbar1', 'hbar2', 'lvsA', 'lvsB', 'lvsC', 'lvsD'] and system_switch['ccsd(t)'] and find_h == 0):
            PR += [xx]
            find_h = 1
        else:
            UCB += [xx]
    if(len(PR) != 1):
        print "Problem1 in independent_mat()!!n\n\n\n"
        printterm(term1)
        pdb.set_trace()
    elif(len(UCB) == 0):
            return [PR, []]
    else:
        flag = 1
        while(flag and len(UCB) > 0):
            AA = independent_mat_help(PR, UCB)
            if(len(AA[0]) == 0):
                flag = 0
            else:
                PR += AA[0]
                UCB = AA[1]
        return [PR, UCB]




def select_up_2_n_B_D(terms, n):       # select terms with 0, 1, 2, .., n-body density matrix elements
    result = []
    for xx in terms:
        flag = 1
        for yy in xx.coefficient.matElement:
            if(yy.name in['nd', 'd1', 'd2', 'd3', 'd4'] and len(yy.matUpperIndicees) > n):
                flag = 0
                break
        if(flag == 1):
            result += [xx]
    return result


def remove_core_cumu(terms):           # remove terms containg cumulants which have core lables
    result = []
    for xx in terms:
        ok = 1
        for yy in xx.coefficient.matElement:
            if(yy.name in ['Lambda', 'lambda', 'lambda1', 'Lambda1']):
                for zz in yy.matUpperIndicees + yy.matLowerIndicees:
                    if(zz.type == 'h' and zz.active == -1):
                        ok = 0
                        break
        if(ok):
            result += [xx]   
    return result




def remove_p_cumu_and_dm(terms):       # remove terms containg cumulants or DMs which have particle lables
    result = []
    for xx in terms:
        ok = 1
        for yy in xx.coefficient.matElement:
            if(yy.name in ['Lambda', 'lambda', 'lambda1', 'Lambda1', 'DEN', 'Gamma', 'gamma', 'd', 'd1', 'd2', 'd3', 'd4']):
                for zz in yy.matUpperIndicees + yy.matLowerIndicees:
                    if(zz.type in ['p', 'e']):
                        ok = 0
                        break
        if(ok): 
            result += [xx]
    return result









#**************************************************
#              comparing & reordering             *
#**************************************************





def fixorder(term1, matname, offset): 

    """ ....A^pq_rs.... -> ....A^qp_sr...., if the two indices spcified by 'offset' don't satisfy (dic_type_value[indA.type] < dic_type_value[indB.type]) """
    """ or (dic_type_value[indA.type] == dic_type_value[indB.type] & int(indA.num) < int(indB.num)).  offset = [0, 1] or offset == [2, 3] """

    term = copy.deepcopy(term1)         
    n = 0                              
    flag = 0
    for i in range(len(term.coefficient.matElement)):
        if(term.coefficient.matElement[i].name == matname):
            j = i
            n += 1
            flag = 1
    if(flag == 0):
        return term1
    else:
        if(n > 1):
            print "\n\n                P1@fixorder.real!"
            printterm(term1)
            print matname
            print offset
            pdb.set_trace()
        else:
            matt = term.coefficient.matElement[j]
            inds = matt.matUpperIndicees + matt.matLowerIndicees
            indA = inds[offset[0]]
            indB = inds[offset[1]]
            if(dic_type_value[indA.type] < dic_type_value[indB.type]):
                return term
            elif(dic_type_value[indA.type] > dic_type_value[indB.type]):
                matt.matUpperIndicees = [matt.matUpperIndicees[1], matt.matUpperIndicees[0]]
                matt.matLowerIndicees = [matt.matLowerIndicees[1], matt.matLowerIndicees[0]]
                return dress_one(term)
            else:
                if(int(indA.num) < int(indB.num)):
                    return term
                elif(int(indA.num) > int(indB.num)):
                    matt.matUpperIndicees = [matt.matUpperIndicees[1], matt.matUpperIndicees[0]]
                    matt.matLowerIndicees = [matt.matLowerIndicees[1], matt.matLowerIndicees[0]]
                    return dress_one(term)
                else:
                    print "\n\n        indA == indB@fixorder.real!"
                    printterm(term1)
                    print matname
                    print offset    
                    pdb.set_trace()


def fixorder_s(terms, matname,offset):
    res = []
    for xx in terms:
        res += [fixorder(xx, matname, offset)]
    return res
                
    

def term_prefer_mats(term, matnames):  # rearrange matelement: put those with names in 'matname' ahead of others. e.g. matnames = ['F_', 't1']
    coll = []
    for i in range(len(matnames)+1):
        coll += [[]] 
    others = []
    for xx in term.coefficient.matElement:
        if(xx.name in matnames):
            coll[matnames.index(xx.name)] += [xx]
        else:
            others += [xx]
    mats = []
    for yy in coll:
        mats += yy
    mats += others
    return Term(Coefficient(term.coefficient.const, mats), term.uOperators)
    



def terms_prefer_mats(terms, matnames): 
    res = []
    for xx in terms:
        res += [term_prefer_mats(xx, matnames)]
    return res



def finalorder(terms):                 # arranged terms according increasing rank, discard terms with  higher than 5
      order0 = []
      order1 = []
      order2 = []
      order3 = []
      order4 = []
      for a in terms:
            if(len(a.uOperators[0].upperIndicees) == 0): order0 = order0 + [a]
            elif(len(a.uOperators[0].upperIndicees) == 1): order1 = order1 + [a]
            elif(len(a.uOperators[0].upperIndicees) == 2): order2 = order2 + [a]
            elif(len(a.uOperators[0].upperIndicees) == 3): order3 = order3 + [a]
            elif(len(a.uOperators[0].upperIndicees) == 4): order4 = order4 + [a]
      return order0 + order1 + order2 + order3+ order4 
      


def dif(list1, list2):                 # returning terms in list2 which do not appear in list1
        res = []
        for x in list2:
                if(x.inArray(list1)):
                    pass
                else:
                    res += [x]
        return res


def indlisteq(list1, list2):           # decide whether the two indlists are equal or not
    res = 1
    if(len(list1) != len(list2)): res = 0
    elif(len(dif(list1, list2)) > 0): res = 0
    elif(len(dif(list2, list1)) > 0): res = 0
    else: res = 1
    return res




def gvhighestpos(numlist):
    res = 0
    num = numlist[0]
    for j in range(len(numlist)):
        if(numlist[j] > num): 
            res = j
    return res



def gvSmallest(numlist):              

    """ return the minimum of values in numlist """ 
    res = numlist[0]
    for i in numlist:
        if(i < res):
            res = i
    return res



def give_highest_indnum(indlist):
    pmax = 0
    emax = 0
    hmax = 0
    mmax = 0
    gmax = 0
    for xx in indlist:
        if(xx.type == 'p' and int(xx.num) > pmax):
            pmax = int(xx.num)
        elif(xx.type == 'e' and int(xx.num) > emax):
            emax = int(xx.num)
        elif(xx.type == 'h' and int(xx.num) > hmax):
            hmax = int(xx.num)
        elif(xx.type == 'm' and int(xx.num) > mmax):
            mmax = int(xx.num)
        elif(xx.type == 'g' and int(xx.num) > gmax):
            gmax = int(xx.num)
    return [pmax, emax, hmax, mmax, gmax]



def rank(list):                        # give the positions of elements of a list in decreasing order
    Num = list
    if(len(list) == 0): 
        return []
    else:
        order = [0]
        for i in range(len(Num)):
            if(i == 0): 
                pass
            else:
                for j in range(len(order)):
                    if(Num[i] <= Num[order[j]]):
                        if(j == len(order)-1):
                            order = order + [i]
                            break
                        else: 
                            pass
                    elif(j == 0):
                        order = [i] + order
                        break
                    else:
                        order = order[0:j] + [i] + order[j:]
                        break
        return order





def rank_indices_by_num(inds):         # order indices of the same type in increasing order of ind.num
    if(len(inds) <= 1):
        return inds
    else:
        indnums = []
        for xx in inds:
            indnums += [int(xx.num)]
        pos = rank(indnums)
        result += []
        for yy in pos:
            result = [inds[pos]] + result
        return result





def pick_1_name(terms, name, pos):     # select those terms: terms[xx].coefficient.matElement[pos].name == name, if not, return terms
    if(len(terms) == 1 or 0):
        return terms
    else:
        result = []
        for ee in terms:
            if(ee.coefficient.matElement[pos].name == name):
                result += [ee]
        if(len(result) > 0):
            return result
        else:
            return terms



def reorder_inds_hmpe(inds):           # reorder indices according to their type h->m->p->e
    hs = []
    Hs = []
    ms = []
    Ms = []
    ps = []
    Ps = []
    es = []
    Es = []
    for xx in inds:
        if(xx.type == 'h'):
            hs += [xx]
        elif(xx.type == 'H'):
            Hs += [xx]
        elif(xx.type == 'm'):
            ms += [xx]
        elif(xx.type == 'M'):
            Ms += [xx]
        elif(xx.type == 'p'):
            ps += [xx]
        elif(xx.type == 'P'):
            Ps += [xx]
        elif(xx.type == 'e'):
            es += [xx]
        elif(xx.type == 'E'):
            Es += [xx]
        else:
            print "\n\n                                   Type unassigned in reorder_inds_hmpe@real.py\n\n"
            pdb.set_trace()
    return hs + Hs + ms + Ms + ps + Ps + es + Es 
    

   

def reorder(term1):                    # arrange indexpairs in operator in an convenient order, pp-pe-...-gg
    if(len(term1.uOperators[0].upperIndicees) < 2):
        return term1
    else:
        mattotal = []
        term = copy.deepcopy(term1)
        for i in range(25):
            mattotal = mattotal + [[]]
        map = []
        totalnum = []
        typelist = ["p", "e", "h", "m", "g"]
        for dd in typelist:
            for ff in typelist:
                map = map + [[dd, ff]]
        oo = term.uOperators[0]
        for ii in range(len(oo.upperIndicees)):
            jj = [oo.upperIndicees[ii].type, oo.lowerIndicees[ii].type]
            mattotal[map.index(jj)] = mattotal[map.index(jj)] + [ii]
        upp = []
        low = []
        for kk in mattotal:
            for ll in kk:
                upp = upp + [oo.upperIndicees[ll]]
                low = low + [oo.lowerIndicees[ll]]
        oo.upperIndicees = upp
        oo.lowerIndicees = low
        term.uOperators[0] = oo
        return term


def reorder_s(terms):
    result = []
    for xx in terms:
        result += [reorder(xx)]
    return result


def max_match(upp_list, low_list):

    """ reorder low_list to maximize matching upp_list and low_list; for non-matching indices, we don't permute indices in low_list """
    """ return reordered low_list """
    result = []
    for yy in range(len(low_list)):
        result += [-1.5]                                   # first fill 'result' list with special numbers         
    help = []
    for xx in low_list:
        if(xx in upp_list):
            result[ upp_list.index(xx) ] = xx              # matching indices 
        else:
            help += [xx]                                   # nonmatching indices
    for xxx in range(len(help)):
        for yyy in range(len(result))[xxx:]:
            if(result[yyy] == -1.5): 
                result[yyy] = help[xxx]                    # fill nonmatching indices one by one in undertemined positions
                break
    return result










#**************************************************
#                   changename                    *
#**************************************************




def terms_changemat_s_name(terms, names1, names2):
    res = []
    for xx in terms: res.append(changemat_s_name(xx, names1, names2))
    return res



def changemat_s_name(term1, names1, names2):
    term = copy.deepcopy(term1)
    namemap = {}
    r = range(len(names1))
    for n in r: namemap[names1[n]] = names2[n] 
    for xx in names1: term = changematname(term, xx, namemap[xx])  # this is different from doing substitution simultenaously
    return term



def changematname(term1, name1, name2):
    term = copy.deepcopy(term1)
    for xx in term.coefficient.matElement:
        if xx.name == name1: xx.name = name2
    return term


def changematname_s(terms, name1, name2):
    result = []
    for xx in terms:
        result += [changematname(xx, name1, name2)]
    return result



def changename(term1):                

    """ change the name of matrixelements, this function was used in x2c dat file, to differentiate different excitation """
    """ operator difference rising from different contractions with Hamiltonian term, refer to x2c file """

    term = copy.deepcopy(term1)    
    if(len(term.coefficient.matElement) == 1):
        return term
    elif(len(term.coefficient.matElement) == 0 or len(term.coefficient.matElement) >  2):
        print "Bug in changename()!!"
        return
    else:
        mat = term.coefficient.matElement[0]
        delta = term.coefficient.matElement[1]
        inau = delta.matUpperIndicees
        inal = delta.matLowerIndicees
        inn = mat.matUpperIndicees + mat.matLowerIndicees
        L = []
        for i in range(len(inau)):
            if(inau[i].inArray(inn)):
                L += [i]
            if(inal[i].inArray(inn)):
                L += [i]
        L = purify(L)
        if(len(inau) == 2):
            if(L == [0]):
                term.coefficient.matElement[1].name = '\Lambda'
            elif(L == [0, 1]):
                pass
            else:
                print "Bug in changename()!!!"
                return
        elif(len(inau) == 3):
            if(L == [0]):
                term.coefficient.matElement[1].name = '\Gamma'
            elif(L == [0, 1]):
                term.coefficient.matElement[1].name = '\Sigma'
            elif(L == [0, 1, 2]):
                pass
            else:
                print "Bug in changename()!!!!"
                return
        elif(len(inau) == 4):
            if(L == [0, 1]):
                term.coefficient.matElement[1].name = '\Omega'
            elif(L == [0, 1, 2]):
                term.coefficient.matElement[1].name = '\Theta'
            elif(L == [0, 1, 2, 3]):
                pass
            else:
                print "Bug in changename()!!!!"
                return



def changename_s(terms):
    result = []
    for xx in terms:
        result += [changename(xx)]
    return result





def changeindname(term1, name1, name2):
    term = copy.deepcopy(term1)
    if(len(term.uOperators)):
        for x in term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees:
            if(x.type == name1):
                x.type = name2
    for xx in term.coefficient.matElement:
        for yy in xx.matUpperIndicees + xx.matLowerIndicees:
            if(yy.type == name1):
                yy.type = name2
    return term


def changeindname_s(terms, name1, name2):
    result = []
    for xx in terms:
        result += [changeindname(xx, name1, name2)]
    return result













#--------------------------------------------------------------------------------------------------------------------------------------

#                                                             Spatial Orbital Canonicalization: connected

#--------------------------------------------------------------------------------------------------------------------------------------









class SF_Connected_Merge_Separator:
    pass









def simplify(terms):
    """ simplifies the array of term objects; adds the terms which are same after canonicalization """

    simplifiedTerms = []
    res = []
    for a in terms:
        inarray = a.inArray2(simplifiedTerms)
        if(inarray[0]):
            termNo = inarray[1]
            simplifiedTerms[termNo].coefficient.const += a.coefficient.const
        else:
            simplifiedTerms.append(a) 
    for b in simplifiedTerms:
        if(abs(b.coefficient.const) > 1e-4): res.append(b)
    return res




#********************************
#              dress            *
#********************************





def makeIndTerm2(term1):               # make ind from a term, first matrixelements then operator
    ind = []
    term = copy.deepcopy(term1)
    for b in term.coefficient.matElement:
        ind = ind + [makeIndexSet(b.matUpperIndicees, b.matLowerIndicees)]
    for a in term.uOperators:
        ind = ind + [makeIndexSet(a.upperIndicees, a.lowerIndicees)]
    return Indicees(ind)


def arrangeIndicees2(ind, term):       # renumber term with ind, begginning with matrixelements, in contrast to arrangeIndicees() 
    for n in range(len(ind.indexSets)):
        set = ind.indexSets[n]
        if(n < len(term.coefficient.matElement)):
            for a in range(len(set.indexPairs)):
                term.coefficient.matElement[n].matUpperIndicees[a] = set.indexPairs[a].upperIndex
                term.coefficient.matElement[n].matLowerIndicees[a] = set.indexPairs[a].lowerIndex
        else:
            for a in range(len(set.indexPairs)):
                term.uOperators[0].upperIndicees[a] = set.indexPairs[a].upperIndex
                term.uOperators[0].lowerIndicees[a] = set.indexPairs[a].lowerIndex
    return term




def dress_one(term):                   # renumber indices of term, in the order of from matrixelements to operator
    x = copy.deepcopy(term)
    ind = makeIndTerm2(x)
    ind = ind.renumber()
    x = arrangeIndicees2(ind, x)
    return x



def dress_all(terms):
    res = []
    for x in terms:
        res += [dress_one(x)]
    return res





#**************************
#      canonical_all      *
#**************************





def term_value(term):                  # calculate the value of a term, we wish it wil give unique value for different orders
    switch = 0
    if(switch):
        pdb.set_trace
    value = 0
    for yy in range(len(term.coefficient.matElement)):
        y = term.coefficient.matElement[yy]
        did = y.matUpperIndicees + y.matLowerIndicees
        for zz in range(len(did)):
            z = did[zz]
            value += (yy + 3) * (zz + 11) * (int(z.num) + 7) * dic_matname_value[y.name] * dic_type_value[z.type] * dic_active_value[z.active] * dic_fix_value[z.fix]
    return value
    

def term_comp_sup(term1, term2):       # in difficult case, we use this select the final one
    res = 0
    flag = 0 
    for xx in range(len(term1.coefficient.matElement)):
        name1 = term1.coefficient.matElement[xx].name
        name2 = term2.coefficient.matElement[xx].name
        if(dic_matname_value[name1] < dic_matname_value[name2]):
            res = 0
            break
        elif(dic_matname_value[name1] > dic_matname_value[name2]):
            res = 1
            break
        else:
            AA = term1.coefficient.matElement[xx].matUpperIndicees + term1.coefficient.matElement[xx].matLowerIndicees
            BB = term2.coefficient.matElement[xx].matUpperIndicees + term2.coefficient.matElement[xx].matLowerIndicees
            for yy in range(len(AA)):
                aa = AA[yy]
                bb = BB[yy]
                if(dic_type_value[aa.type] < dic_type_value[bb.type]):
                    res = 0
                    flag = 1
                    break
                elif(dic_type_value[aa.type] > dic_type_value[bb.type]):
                    res = 1
                    flag = 1
                    break
                elif(int(aa.num) < int(bb.num)):
                    res = 0
                    flag = 1
                    break
                elif(int(aa.num) > int(bb.num)):
                    res = 1
                    flag = 1
                    break
        if(flag):
            break
    return res

            
    


def term_value_2(term):                # calculate the value of a term, we wish it wil give unique value for different orders, NOT IN USE
    switch = 0
    if(switch):
        pdb.set_trace
    value = 0
    overlap = []
    for yy in range(len(term.coefficient.matElement)):
        y = term.coefficient.matElement[yy]
        did = y.matUpperIndicees + y.matLowerIndicees
        for zz in range(len(did)):
            z = did[zz]
            if(z.inArray(overlap)):
                value += (yy + 3) * (zz + 11) * (int(z.num) + 7) * dic_matname_value[y.name] * dic_type_value[z.type] * dic_active_value[z.active] * dic_fix_value[z.fix]
            else:
                value += (yy + 3) * (zz + 11) * (int(z.num) + 219) * dic_matname_value[y.name] * dic_type_value[z.type] * dic_active_value[z.active] * dic_fix_value[z.fix]
                overlap += [z]
    return value




def mat_putfirst(mat, pos):            # decide the order of the 2(or 1) indexpairs, it works for rank2, so be careful when rank >= 3 
    if(len(mat.matUpperIndicees) > 2):
        print "\n\n                    Not implemented for three or higher-body operators in mat_putfirst@real.py\n\n" 
        pdb.set_trace()
        return
    elif(mat.name in ['s3', 's4', 'w']):
        return mat
    elif(len(mat.matUpperIndicees) == 1):
        return mat
    elif(pos == 0 or pos == 3):
        #print "\n\n pos = %s in mat_putfirst() \n\n" % pos
        return mat
    else:
        return MatElement(mat.name, [mat.matUpperIndicees[1], mat.matUpperIndicees[0]], [mat.matLowerIndicees[1], mat.matLowerIndicees[0]], mat.connectivity, mat.other) 
    


def canonical_one_help(mats1, mats2):  

    """ divide mats2 to 2 classes according to mats1, and sort them, by the order of appering in indices of mats1, """
    """ return B1 and B2; Looks like this is the only place to change the order of index pairs """
 
    NY = []                      
    NY3 = []
    NY6 = []
    NY7 = []
    WS = []
    if(len(mats2) == 0):
        return [[], []] 
    else:
        ind = []
        for x in mats1:
            ind += x.matUpperIndicees + x.matLowerIndicees
        for xx in mats2:
            NY2 = []
            NY5 = []
            flag = 0
            ind2 = xx.matUpperIndicees + xx.matLowerIndicees
            for yy in range(len(ind2)):
                if(ind2[yy].inArray(ind)):
                    flag = 1
                    NY5 += [yy]
                    NY2 += [ind2[yy].gvIndexNo(ind)]
            if(flag):
                NW = rank(NY2)[-1]
                NY3 += [mat_putfirst(xx, NY5[NW])] # due to the possible confusion in 'w' when getting its indices' original order,  
                NY7 += [NY2[NW]]                   # we handle it specially. Check mat_putfirst()!!
            else:
                WS += [xx]
        NY4 = rank(NY7)
        for No in range(len(NY4)):
            NK = (-1)*(No + 1)
            NY += [NY3[NY4[NK]]] 
        return [NY, WS]  



def canonical_one(term1):
    term = copy.deepcopy(term1)
    PR = []
    UCB = []
    find_h = 0
    for xx in term.coefficient.matElement:
        if(xx.name in ['f', 'v', 'nf', 'nv', 'kmf', 'kmv', 'hbar1', 'hbar2', 'lvsA', 'lvsB', 'lvsC', 'lvsD'] and system_switch['ccsd(t)'] == 0):
            PR += [xx]
        elif(xx.name in ['f', 'v', 'nf', 'nv', 'kmf', 'kmv', 'hbar1', 'hbar2', 'lvsA', 'lvsB', 'lvsC', 'lvsD'] and system_switch['ccsd(t)'] and find_h == 0):
            PR += [xx]                 # we handle ccsd(t) separately because 'v' appear more than once  
            find_h = 1
        else:
            UCB += [xx]
    if(len(PR) != 1):
        print "Problem1 in canonical_one()!!n\n\n\n"
        printterm(term)
        pdb.set_trace()
    elif(len(UCB) == 0 and PR[0].name in ['f', 'nf', 'kmf', 'hbar1']):
        return dress_one(term)     # we do 'dress' here!!!!!!!!!!!!!!!!!!!!!!!!!!
    else:
        if(PR[0].name in ['f', 'nf', 'kmf', 'hbar1']):
            termm = Term(Coefficient(term.coefficient.const, PR), term.uOperators)          
            flag = 1
            while(flag):
                AA = canonical_one_help(PR, UCB)
                termm.coefficient.matElement += AA[0]
                PR = copy.deepcopy(AA[0])
                UCB = copy.deepcopy(AA[1])
                if(len(UCB) == 0):
                    flag = 0
            return dress_one(termm)    
        else:
            res = []
            mat1 = copy.deepcopy(PR[0])
            mat2 = MatElement(mat1.name, [mat1.matUpperIndicees[1], mat1.matUpperIndicees[0]], [mat1.matLowerIndicees[1], mat1.matLowerIndicees[0]], mat1.connectivity, mat1.other)
            vs = [[mat1], [mat2]]
            for i in range(2):
                IAS = copy.deepcopy(vs[i])
                PEN = copy.deepcopy(UCB)
                termm = Term(Coefficient(term.coefficient.const, IAS), term.uOperators, term.type, term.other)
                flag = 1
                while(flag):
                    AA = canonical_one_help(IAS, PEN)
                    termm.coefficient.matElement += AA[0]
                    IAS = copy.deepcopy(AA[0])
                    PEN = copy.deepcopy(AA[1])
                    if(len(PEN) == 0):                       
                        flag = 0
                    elif(len(IAS) == 0):
                        print "Problem2 in canonical_one()!\n\n"               # Some matrixelement has no contraction with the others
                        pdb.set_trace()                            
                res += [dress_one(termm)]
            if(term_equal(res[0], res[1])):
                return res[0]         
            else:                      # Here  we deal with the permuted form of 'v', in two steps  
                switch1 = 1            
                if(switch1):
                    aa = term_value(res[0])
                    bb = term_value(res[1])
                    if(aa < bb):
                        return res[0]
                    elif(aa > bb):
                        return res[1]
                    else:
                        #print "term value equal \n\n"
                        switch2 = 0
                        if(switch2):
                            pdb.set_trace()
                            return res[0]
                switch3 = 1
                if(switch3):
                    n = term_comp_sup(res[0], res[1])
                    return res[n]



def canonical_all(terms):              # if merge_all() gets stuck, then check canonical_all!!! Probably there are uncontracted mat in some terms, then
    res = []                           # activate system_switch['MCSCF_mi'] to deal with it
    I = 1
    printonce = 1
    for xx in terms:
        if(printonce):
            print "\n              Take care of possible uncontracted mat\n"
            printonce = 0
        inde = independent_mat(xx)
        if(len(inde[1]) == 0):
            res += [canonical_one(xx)]
        else:
            print "\n     disconnected terms appear, change the canonical scheme to 'p' \n"
            printterm(xx)
            pdb.set_trace()
        I += 1
        if(I - 100 * (I/100) == 0 or 0):
            print "\n\n                                                     I = %s" % I
    return res







#**************************
#        merge_all        *
#**************************





def term_equal(term1, term2):
    flag = 1
    if(term1.uOperators[0].equals(term2.uOperators[0]) == 0):
        flag = 0
    if(flag):
        if(len(term1.coefficient.matElement) != len(term2.coefficient.matElement)):
            flag = 0
        else:
            for ss in range(len(term1.coefficient.matElement)):
                if(term1.coefficient.matElement[ss].name != term2.coefficient.matElement[ss].name):
                    flag = 0
                    break
                elif(len(term1.coefficient.matElement[ss].matUpperIndicees) != len(term2.coefficient.matElement[ss].matLowerIndicees)):
                    flag = 0
                    break
                else:
                    for qq in range(len(term1.coefficient.matElement[ss].matUpperIndicees)):
                        if(term1.coefficient.matElement[ss].matUpperIndicees[qq].equal(term2.coefficient.matElement[ss].matUpperIndicees[qq]) == 0):
                            flag = 0
                            break
                        elif(term1.coefficient.matElement[ss].matLowerIndicees[qq].equal(term2.coefficient.matElement[ss].matLowerIndicees[qq]) == 0):
                            flag = 0
                            break
                    if(flag == 0):
                        break
    return flag


def delete_mat(term, matnames):         # delete matrixelement of name in matnames from 'term'
    mats = []
    for xx in term.coefficient.matElement:
        if(xx.name in matnames):
            pass
        else:
            mats += [xx]
    return Term(Coefficient(term.coefficient.const, mats), term.uOperators)



def merge_all(terms1):                                              # if merge_all() gets stuck, then check canonical_all!!! Probably there are uncontracted mat       
    if(system_switch['use spin orbital'] == 1):                     # in some terms, then activate system_switch['MCSCF_mi'] to deal with it
        return so_merge_all(terms1)
    elif(system_switch['canonical scheme'] == 'p'):
        return merge_all_by_permu(terms1)
    elif(system_switch['canonical scheme'] != 'c'):
        print "\n undefined canonical scheme in merge_all() \n"
        pdb.set_trace()
    else:
        terms = copy.deepcopy(terms1)
        terms = canonical_all(terms)                                                   
        print "                                             canonical_all() done. \n"   # Here we first canonicalize all terms, then merge them
        print "                                             simplication starts: " + str(len(terms)) + " terms. \n"
        res = []
        result = []
        counter = 0
        for xx in terms:
            counter += 1
            if(counter - 200 * (counter/200) == 0): print "                                             counter = %s \n" % counter
            if(len(result) == 0): result += [xx]
            else:
                flag = 0
                for yy in result:
                    if(term_equal(delete_mat(xx, ['X']), delete_mat(yy, ['X']))):  # because 'X' matrixelement is virtual, so it should not bother the comparison
                        yy.coefficient.const += xx.coefficient.const
                        flag = 1
                        break
                if(flag == 0): result += [xx]
        for zz in result:
            if(abs(zz.coefficient.const) > 0.00001): res.append(zz) 
        print "                                             simplication done: " + str(len(res)) + " terms. \n"
        return res





def merge_all_2_parts(terms1, terms2):# assume that terms1 and terms2 are already canonicalized, we simplify terms1 + terms2
    terms = copy.deepcopy(terms2)        
    print "                                             simplication starts.  \n"
    print "                                             No. of terms = %s\n" % len(terms)
    res = []
    result = copy.deepcopy(terms1)
    counter = 0
    for xx in terms:
        counter += 1
        if(counter - 100 * (counter/100) == 0):
            print "                                             counter = %s \n" % counter
        if(len(result) == 0):
            result += [xx]
        else:
            flag = 0
            for yy in result:
                if(term_equal(delete_mat(xx, ['X']), delete_mat(yy, ['X']))):  # because 'X' matrixelement is virtual, so it should not bother the comparison
                    yy.coefficient.const += xx.coefficient.const
                    flag = 1
                    break
            if(flag == 0):
                result += [xx]
    for zz in result:
        if(abs(zz.coefficient.const) > 0.00001):
            res += [zz]
    print "                                             simplication done: " + str(len(res)) + " terms. \n"

    return res



















#--------------------------------------------------------------------------------------------------------------------------------------

#                                                             Spatial Orbital Canonicalization: permutation

#--------------------------------------------------------------------------------------------------------------------------------------










class SF_Permutation_Merge_Separator:
    def __init__(self, start = "Y"):
        self.start = type















#***********************************************
#   Find Possible Permuted Forms of One Term   *
#***********************************************






#  permute mats -> permute index pairs of each mat 






def mat_2_mats(mat):                        # given a 'mat', we read its indlist, then permute, screen, contruct possible mats
    indset = mat.mat_2_indexset()
    indpair_s_permu = permu_general(indset.indexPairs)
    indlists = []
    for AA in indpair_s_permu:
        indlists += [IndexSet(AA).get_list()]
    return screen_indlists_to_mats(indlists, mat.name)



def mats_2_matslist(mats):                  # given 'mats' (from term.coefficient.matElement), return possible mats by doing permutation and screening
    A = []                                  
    for xx in mats:
        A += [mat_2_mats(xx)]
    return combi_lists(A)



def term_from_permu_2_terms_help(term):         
    result = []                           
    A = mats_2_matslist(term.coefficient.matElement)   
    if(0):
        print "\n     Step 1 in term_from_permu_2_terms done\n" 
    const = term.coefficient.const
    for xx in A:
        result += [Term(Coefficient(const, xx), term.uOperators)]
    return result
   

def term_from_permu_2_terms(term):
    """ given a term, return possible terms by doing permutation of index pairs in term.coefficient.matElement and screening, return terms """

    if len(term.coefficient.matElement) == 0: return [term]
    else: return term_from_permu_2_terms_help(term)
 


def term_variants(term):
    #set_trace()
    terms = term_permu_mat(term)                                               # permute matelemenmts; in this new version, the mats are automatically grouped and ordered according to   
    A = []                                                                     # dic_matname_value 
    for xx in terms: A += term_from_permu_2_terms(xx)                          # permute indices, '2' means 'to'
    return A





#*************************************************
#   Combine Everyting to One Canonical Function  *
#*************************************************





def canonical_one_by_permu(term):                                              # first permute mats and their index pairs, then select the one with the highest priority
    terms = term_variants(term)                                                # the underlying assumption is that (1) all variants have the the same operator part (indeed)
    result = compare_all_term_variants(terms)                                  # (2) so keep operator part unchanged; permute matElement part and then identify the unique
    result1 = dress_one(result)                                                #     matElement. Finally do 'dress_one' to the matElement and operator !!! 
    return result1


def canonical_all_by_permu(terms):
    print "\n                                             canonical_all_by_permu starts\n"
    result = []
    clock = 1
    checknum = 100
    for xx in terms:
        result.append(canonical_one_by_permu(xx))
        if(clock - checknum *(clock/checknum) == 0 and 1): print "                                             processing term %s \n"  % clock
        clock += 1
    print "                                             canonical_all_by_permu done. \n"   
    return result




    
def merge_all_by_permu(terms1):            
    terms = copy.deepcopy(terms1)      
    terms = canonical_all_by_permu(terms) # Here we first canonicalize all terms, then merge them
    print "                                             simplication starts: " + str(len(terms)) + " terms. \n"
    res = []
    result = []
    counter = 0
    for xx in terms:
        counter += 1
        if(counter - 100 * (counter/100) == 0): print "                                             processing term %s \n" % counter
        if(len(result) == 0): result += [xx]
        else:
            flag = 0
            for yy in result:
                if(term_equal(delete_mat(xx, ['X']), delete_mat(yy, ['X']))):  # because 'X' matrixelement is virtual, so it should not bother the comparison
                    yy.coefficient.const += xx.coefficient.const
                    flag = 1
                    break
            if(flag == 0): result += [xx]
    for zz in result:
        if(abs(zz.coefficient.const) > 0.00001): res.append(zz) 
    print "                                             simplification done: " + str(len(res)) + " terms. \n"
    return res
























#--------------------------------------------------------------------------------------------------------------------------------------

#                                                    Spin Orbital Canonicalization

#--------------------------------------------------------------------------------------------------------------------------------------











class SoMergeSeparator:
    def __init__(self, start = "Y"):
        self.start = type













#***********************************************
#                    Parity                    *
#***********************************************





""" Take care of the sign change by monitoring the change of the order of indices in op and mats """



def list_parity_change(Indlist1, Indlist2): 

    """ determine the change of the sign due to change of the order of the indices """

    indlist1 = copy.deepcopy(Indlist1)
    indlist2 = copy.deepcopy(Indlist2)
    result = 1                              # return 1 or -1
    for xx in range(len(indlist1)):
        indA = indlist1[xx]
        if(indA.inArray(indlist2) == 0):
            print "\n\n                     The ind unfound in indlist2. list_parity_change@real.py\n\n"
            print indA.gvIndex()
            print "\n\n                     The first index list:\n"
            showindlist(indlist1)
            print "\n\n                     The second index list:\n"
            showindlist(indlist2)  
            pdb.set_trace()
        if(indA.equal(indlist2[xx])):
            pass
        else:
            pos = indA.gvIndexNo(indlist2)
            indB = indlist2[xx]
            indlist2[xx] = indA
            indlist2[pos] = indB
            result *= -1
    return result 



def mat_parity_change(mat1, mat2):          
    
    """ determine the change of the sign due to change of the order of upper and lower indices, return 1 or -1 """
  
    result = 1                        
    result *= list_parity_change(mat1.matUpperIndicees, mat2.matUpperIndicees)
    result *= list_parity_change(mat1.matLowerIndicees, mat2.matLowerIndicees)
    return result



def mats_parity_change(mats1, mats2):
    result = 1
    for xx in range(len(mats1)):
        result *= mat_parity_change(mats1[xx], mats2[xx])
    return result



def term_parity_change_help_find_order(term1, term2):                

    """ after permutation, the order of mats changed; to facilitate doing term_parity_change(), we make them in the same order """        

    mats_inter = []                                                   
    for xx in term2.coefficient.matElement:
        flag = 0
        for yy in term1.coefficient.matElement:
            if(yy.name == xx.name and len(dif(yy.matUpperIndicees, xx.matUpperIndicees)) == 0 and len(dif(yy.matLowerIndicees, xx.matLowerIndicees)) == 0):
                flag = 1
                mats_inter += [yy]
                break
        if(flag == 0):
            print "\n\n                     Not match? term_parity_change_help_find_order@real.py\n\n"
            printterms([term1, term2])
            pdb.set_trace()
    return Term(Coefficient(term1.coefficient.const, mats_inter), term1.uOperators)  



def term_parity_change(termA, term2):       # term1 is the original one
    term1 = term_parity_change_help_find_order(termA, term2)
    result = 1
    result *= mats_parity_change(term1.coefficient.matElement, term2.coefficient.matElement)
    result *= list_parity_change(term1.uOperators[0].upperIndicees, term2.uOperators[0].upperIndicees)
    result *= list_parity_change(term1.uOperators[0].lowerIndicees, term2.uOperators[0].lowerIndicees)
    return result
   
 



#***********************************************
#   Find Possible Permuted Forms of One Term   *
#***********************************************






#  permute mats -> permute indices of each mat -> for each index-pemuted mat, give the right order of indices in operator


def select_ind_from_2inds(ind1, ind2):        
    if(1):
        if ind1.type != ind2.type: 
            if ind1.type in ['P', 'E', 'H', 'M'] and ind2.type in ['p', 'e', 'h', 'm']: return 0       # put alpha before beta indices, to allow only alpha-beta type of quantities
            elif ind2.type in ['P', 'E', 'H', 'M'] and ind1.type in ['p', 'e', 'h', 'm']: return 1
            elif dic_type_value[ind1.type] < dic_type_value[ind2.type]: return 0
            elif dic_type_value[ind1.type] > dic_type_value[ind2.type]: return 1
            else: set_trace()
        else:
            if ind1.fix != ind2.fix:
                if ind1.fix > ind2.fix: return 0
                else: return 1
            elif ind1.fix == ind2.fix == 1 and int(ind1.num) < int(ind2.num): return 0
            elif ind1.fix == ind2.fix == 1 and int(ind1.num) > int(ind2.num): return 1
            elif ind1.active != ind2.active :
                if ind1.active < ind2.active: return 0
                elif ind1.active > ind2.active: return 1
                else: set_trace()
            else:
                return 3



def select_indlist_2(inds1, inds2):        
    """ select the inds prior to the other one, if they have the same priority, return 2 """
    """ pay attention to the order of selecting, which should be the same as direct_compare_2_mats """ 

    LL = range(len(inds1))
    for m in LL:
        res = select_ind_from_2inds(inds1[m], inds2[m])
        if res != 3: return res
    return 3



def select_indlists(indlists):             # return the indlist which is prior to the others
    if(len(indlists) == 0):
        print "\n\n                         P1: select_indliists@real.py\n\n"
        pdb.set_trace()
        return
    elif(len(indlists) == 1):
        return indlists
    else:
        result = [indlists[0]]
        for M in range(len(indlists)):
            if(M < len(indlists) -1):
                pos = select_indlist_2(result[-1], indlists[M + 1])
                if(pos == 0):
                    pass
                elif(pos == 1):
                    result = [indlists[M+1]]
                elif(pos == 3):
                    result += [indlists[M+1]]
                else:
                    print "\n\n             P2: select_indliists@real.py\n\n"
                    pdb.set_trace()
                    return
        return result



def indlist_to_mat(indlist, matname):       # construct a matelement from an index list
    return MatElement(matname, indlist[:len(indlist)/2], indlist[len(indlist)/2:])



def indlists_to_mats(indlists, matname):    # refer to the above  
    result = []
    for xx in indlists:
        result += [indlist_to_mat(xx, matname)]
    return result



def screen_indlists_to_mats(indlists, matname):                                # given 'indlists', we screen them and contruct possible mats
    return indlists_to_mats(select_indlists(indlists), matname)



def so_mat_2_mats(mat):                        # given a 'mat', we read its indlist, then permute, screen, contruct possible mats
    inds1 = mat.matUpperIndicees
    inds2 = mat.matLowerIndicees
    p_inds1 = permu_general(inds1)
    p_inds2 = permu_general(inds2)
    indlists = combi_lists_2([p_inds1, p_inds2])
    return screen_indlists_to_mats(indlists, mat.name)



def so_mats_2_matslist(mats):                  # given 'mats' (from term.coefficient.matElement), return possible mats by doing permutation and screening
    A = []
    for xx in mats:
        A += [so_mat_2_mats(xx)]
    return combi_lists(A)



def change_op_acc_mats_help_1(inds1, inds2, indsfrommats): 

    """ indsfrommats are indices from term.coefficient.matElement; 'inds1' and 'inds2' are from permuations of indices of Uoperator, """
    """ compare with indsfrommats, return 0 or 1 corresponding to the one with higher priority """

    LIST = indsfrommats                
    result = 3                        
    for m in range(len(inds1)):
        result = select_indlist_2([inds1[m]], [inds2[m]])
        if(result == 3):
            within1 = inds1[m].inArray(LIST)
            within2 = inds2[m].inArray(LIST)
            if(within1 == 0 and within2 == 1):
                result = 0
                break
            elif(within1 == 1 and within2 == 0):
                result = 1
                break
            elif(within1 == 1 and within2 == 1):
                id1 = inds1[m].gvIndexNo(LIST)
                id2 = inds2[m].gvIndexNo(LIST)
                if(id1 < id2):
                    result = 0
                    break
                elif(id1 > id2):
                    result = 1
                    break
                else:
                    pass
            else:
                pass
        if(result in [0, 1]):
            break
        else:
            pass
    return result



def change_op_acc_mats_help_2(indlists, indsfrommats):                    

    """  given indlists, select the one with the highest priority by comparing with 'indsfrommats """

    if(len(indlists) == 0):   
        print "\n\n                      P1: this case should be avoided in change_op_acc_mats@real.py\n\n"
        pdb.set_trace()
    elif(len(indlists) == 1):
        return indlists[0]
    else:
        result = indlists[0]
        for M in range(len(indlists)):
            if(M < len(indlists) -1):
                pos = change_op_acc_mats_help_1(result, indlists[M + 1], indsfrommats)
                if(pos in [0, 3]):
                    pass
                elif(pos == 1):
                    result = indlists[M+1]
                else:
                    print "\n\n             P2: select_indliists@real.py\n\n"
                    pdb.set_trace()
        return result



def change_op_acc_mats(op, mats):        

    """  adjust the order of indices in operator according to mats, return [Uoperator] """

    indlist = op.upperIndicees + op.lowerIndicees
    indsfrommats = []
    for xx in mats:
        for yy in xx.matUpperIndicees + xx.matLowerIndicees:
            if(yy.inArray(indsfrommats) == 0):
                indsfrommats += [yy]
    indlists = combi_lists_2([permu_general(op.upperIndicees), permu_general(op.lowerIndicees)]) 
    if(len(indlists) > 0):
        opindlist = change_op_acc_mats_help_2(indlists, indsfrommats)
    else:
        opindlist = []
    if(len(opindlist) > 0):
        result = [Uoperator(opindlist[:len(opindlist)/2], opindlist[len(opindlist)/2:])] 
    else:
        result = [Uoperator([], [])]
    return result    



def so_term_from_permu_2_terms(term):         

    """ given a term , return possible terms by doing permutation of indices in term.coefficient.matElement and screening, return terms """

    if(len(term.coefficient.matElement) == 0):
        return [term]
    else:
        result = []                           
        A = so_mats_2_matslist(term.coefficient.matElement)
        if(0):
            print "\n     Step 1 in so_term_from_permu_2_terms done\n" 
        const = term.coefficient.const
        for xx in A:
            result += [Term(Coefficient(const, xx), term.uOperators)]
        return result
   
 
def term_permu_mat(term):
    """ permute term.coefficient.matElement, return terms """

    if len(term.coefficient.matElement): return term_permu_mat_help(term)
    else: return [term]


def term_permu_mat_help(term):                  
    #set_trace()
    result = []
    Vf = []     # it is easiest to assign a list to each individual type of mat, even if some of them are exclusive to each other
    Vv = []
    Vnf = []
    Vnv = []
    Vkmf = []
    Vkmv = []
    Vhbar1 = []
    Vhbar2 = []
    Vg = []
    Vk = []
    T1 = []
    T2 = []
    T3 = []
    S1 = []
    S2 = []
    S3 = []
    S4 = []
    W = []
    Dd1 = []
    Dd2 = []
    Dd3 = []
    Dd4 = []
    X1 = []
    X2 = []
    X3 = []
    ND = []
    X = []
    DELTA = []
    gamma = []
    xambda = []
    xambda1 = []
    Gamma = []
    Lambda = []
    Lambda1 = []
    eta = []
    DEN = []
    Cbar = []
    DI = []
    DbI = []
    r12 = []
    matname_2_container = {'f': Vf, 'v': Vv, 'nf': Vnf, 'nv': Vnv, 'kmf': Vkmf, 'kmv': Vkmv, 'hbar1': Vhbar1, 'hbar2': Vhbar2, 'G': Vg, 'K': Vk, 't1': T1, 't2': T2, 't3': T3, 's1': S1, 's2': S2, 's3': S3, 's4': S4, 'w': W, 'd1': Dd1, 'd2': Dd2, 'd3': Dd3, 'd4': Dd4, 'x1': X1, 'x2': X2, 'x3': X3, 'nd': ND, 'X': X, 'delta': DELTA, 'lambda': xambda, 'lambda1': xambda1, 'Lambda': Lambda, 'Lambda1': Lambda1, 'gamma': gamma, 'Gamma': Gamma, 'eta': eta, 'DEN': DEN, 'Cbar': Cbar, 'D_INV': DI, 'Cbar_INV': DbI, 'r': r12}
    #set_trace()
    for xx in term.coefficient.matElement: matname_2_container[xx.name].append(xx)
    #set_trace()
    alli = dic_matname_value.items()
    items = [(v, k) for k, v in alli]
    items.sort()
    items = [(k, v) for v, k in items]   # this few lines sorts the keys according to the values
    #set_trace()
    A = []
    for yy in items:
        container = matname_2_container[yy[0]]
        if len(container): A.append(permu_general(container))
    #set_trace()
    B = combi_lists_2(A)
    C = []
    for yy in B: C.append(Term(Coefficient(term.coefficient.const, yy), term.uOperators)) 
    #set_trace()
    return C
    


def so_term_variants(term):
    terms = term_permu_mat(term)                                               # permute matelemenmts  
    A = []
    for xx in terms:
        A += so_term_from_permu_2_terms(xx)                                       # permute indices
    return A





#**************************************************
#   Compare Permuted Forms, Find the Unique One   *
#**************************************************





""" The idea here is to compare two terms in this order: """
"""   1 we only compare the mats of term.coefficient.matElement one by one """
"""   2 P(mat with lower dic_mat_value) > P(mat with higher dic_mat_value) """
"""   3 When compare corresponding indices of mats, first upper then lower """
"""   4 P > E > H > M > p > e > h > m """
"""   5 P(fixed index) > P(free index) """
"""   6 P(external index) > P(contracted index) """
"""   7 P(index appearing earlier in contracted index list) > P(index appearing later in contracted index list) """

""" This scheme works because: for fixed indices, we can find a unique order for them; for contracted indices, only their 'types' and positions, not 'num',  matters; """
"""                            for free but not fixed indices, only 'types' matters """
""" This scheme can be modified for spatial canonicalization """



def direct_compare_2_mats(mats1, mats2): 
    """ compare mats1 and mats2, find the one which survives, return 0 or 1 """

    result = 3
    LIST_1 = []
    LIST_2 = []
    flag = 'N'                              # 'N' means 'no'; once comparison ends, 
    if len(mats1) != len(mats2):
       printMatElements(mats1)
       printMatElements(mats2)
       pdb.set_trace()
    for n in range(len(mats1)):
        if dic_matname_value[mats1[n].name] < dic_matname_value[mats2[n].name]:                                       # first compare 'matname'    
            result = 0
            break
        elif dic_matname_value[mats1[n].name] > dic_matname_value[mats2[n].name]:
            result = 1
            break
        else: 
            inds1 = mats1[n].matUpperIndicees + mats1[n].matLowerIndicees
            inds2 = mats2[n].matUpperIndicees + mats2[n].matLowerIndicees
            for m in range(len(inds1)):
                result = select_ind_from_2inds(inds1[m], inds2[m])
                if result == 3:                                                                                       # preliminary comparison can not differentiate
                    within1 = inds1[m].inArray(LIST_1)
                    within2 = inds2[m].inArray(LIST_2)
                    if within1 == 0 and within2 == 1:
                        result = 0
                        break
                    elif within1 == 1 and within2 == 0:
                        result = 1
                        break
                    elif within1 == 1 and within2 == 1:
                        id1 = inds1[m].gvIndexNo(LIST_1)
                        id2 = inds2[m].gvIndexNo(LIST_2)
                        if id1 < id2:
                            result = 0
                            break
                        elif id1 > id2:
                            result = 1
                            break
                        else: pass
                    else:
                        LIST_1 += [inds1[m]]
                        LIST_2 += [inds2[m]]
                        pass
                if result in [0, 1]: break                                                                                  # if direct inds comparison can differentiate the inds
                else: pass         
        if result in [0, 1]: break
        else: pass
    if result in [0, 1]: return result
    else: return 0
    


def direct_compare_2_term_variant(term1, term2):                   
    """ compare term1 and term2, find the one which survives  """

    result = direct_compare_2_mats(term1.coefficient.matElement, term2.coefficient.matElement)
    if result == 0: return term1
    elif result == 1: return term2
    else:
        print result
        set_trace() 
        return



def compare_all_term_variants(terms):
    if(len(terms) == 0): set_trace()
    elif len(terms) == 1: return terms[0]
    else: 
        result = terms[0]
        for xx in range(len(terms)):
            if xx < len(terms)-1:
                result = direct_compare_2_term_variant(result, terms[xx + 1])
        return result





#*************************************************
#   Combine Everyting to One Canonical Function  *
#*************************************************





def so_canonical_one(term):                                                    # first permute mats and their indices, then select  the one with the highest
    terms = so_term_variants(term)                                             # priority, then determine the sign change
    result = compare_all_term_variants(terms)
    result2 = Term(result.coefficient, change_op_acc_mats(result.uOperators[0], result.coefficient.matElement))
    result2.coefficient.const *= term_parity_change(term, result2)
    #if(system_switch['so_canonical_one_do_dress'] == 1):
    result3 = dress_one(result2)
    #else:
    #result3 = result2
    return result3


def so_canonical_all(terms):
    result = []
    clock = 1
    for xx in terms:
        if(clock - 10 *(clock/10) == 0 and 0):
            print "\n   the %s th term in so_canonical_all\n" % clock
        result += [so_canonical_one(xx)]
        clock += 1
    return result



def so_merge_all(terms1):              
    terms = copy.deepcopy(terms1)      
    terms = so_canonical_all(terms)
    print "                                             so_canonical_all() done. \n"   # Here we first canonicalize all terms, then merge them
    print "                                             simplication starts: " + str(len(terms)) + " terms. \n" 
    res = []
    result = []
    counter = 0
    for xx in terms: 
        counter += 1
        if(counter - 200 * (counter/200) == 0): print "                                             counter = %s \n" % counter
        if(len(result) == 0): result += [xx]
        else:
            flag = 0 
            for yy in result:
                if(term_equal(delete_mat(xx, ['X']), delete_mat(yy, ['X']))):  # because 'X' matrixelement is virtual, so it should not bother the comparison
                    yy.coefficient.const += xx.coefficient.const
                    flag = 1
                    break
            if(flag == 0): result += [xx]
    for zz in result:
        if(abs(zz.coefficient.const) > 0.00001): res.append(zz)
    print "                                             simplification done: " + str(len(res)) + " terms. \n"
    return res






























#--------------------------------------------------------------------------------------------------------------------------------------

#                                               spin expansion

#--------------------------------------------------------------------------------------------------------------------------------------











class ExplicitSpinSeparator:
    pass





# assume that operator should not change Sz and there is no 'g' index anywhere




def explicit_spin_terms(terms):   # if a term has different number of alpha(beta) indices in the upper and lower indices of the operator or in
    result = []                   # any matelement, it is discarded from the expansion
    if 1:
        print str(len(terms)) + " terms@explicit_spin_terms\n"
        count = 1
    for xx in terms: 
        result.extend(explicit_spin_term(xx))
        if 1:
            if count % 10 == 0: print str(count) + "-th term,\n"
            count += 1
    return result



def explicit_spin_term(term):     # write one term in explicit spin forms; first mat elements indices, then operator indices. 
    terms = explicit_spin_term_mats(term)
    result = []
    for xx in terms: result += explicit_spin_term_operator(xx)
    return result



def explicit_spin_term_operator(term):
    A = []
    indlist = []
    for xx in term.uOperators[0].upperIndicees + term.uOperators[0].lowerIndicees:
        if(xx.other[0] == ""):
            indlist += [xx]
            A += [[Index(xx.type, xx.num, xx.att, xx.fix, xx.active, ['explicit'] + xx.other[1:]), Index(dic_beta_to_alpha[xx.type], xx.num, xx.att, xx.fix, xx.active, ['explicit'] + xx.other[1:])]]
        elif(xx.other[0] == 'explicit'):
            pass
        else:
            print "\n\n                     index.other[0] unexpected: explicit_spin_term_operator@real.py\n\n"
            pdb.set_trace()
    if(len(indlist) == 0):
        return [term]
    else:
        B = combi_lists(A)
        C = []
        for yy in B:
            C += [reset(term, indlist, yy)]
        D = []
        for zz in C:
            if(explicit_spin_check(zz, 0) == 1):# make sure there are the same number of alpha and beta indices; the second argument is set to zero 
                D += [zz]                       # because up to here there should be no 'general' indices
        return D



def explicit_spin_term_mats(term):           # write term.coefficient.matElement in explicit spin
    result = [term]
    for xx in range(len(term.coefficient.matElement)):
        result = explicit_spin_terms_mat(result, xx)
    return result



def explicit_spin_terms_mat(terms, n):      # write the n-th mat element of terms in explicit spin 
    result = []
    for xx in terms:
        result += explicit_spin_term_mat(xx, n)
    return result



dic_beta_to_alpha = {'p': 'P',
                     'e': 'E',
                     'h': 'H',
                     'm': 'M'}



def explicit_spin_term_mat(term, n):        # write the n-th matelement in explicit spin; assume that in any single mat element, there are no repeated indices
    A = []
    indlist = []
    for xx in term.coefficient.matElement[n].matUpperIndicees + term.coefficient.matElement[n].matLowerIndicees:
        if(xx.other[0] == ""):
            indlist += [xx]
            A += [[Index(xx.type, xx.num, xx.att, xx.fix, xx.active, ['explicit'] + xx.other[1:]), Index(dic_beta_to_alpha[xx.type], xx.num, xx.att, xx.fix, xx.active, ['explicit'] + xx.other[1:])]]
        elif(xx.other[0] == 'explicit'):
            pass
        else:
            print "\n\n                     index.other[0] unexpected: explicit_spin_term_mat@real.py\n\n"
            pdb.set_trace()
    if(len(indlist) == 0):
        return [term]
    else:
        B = combi_lists(A)                 # pick up one element from each A element and put them in a list, which is an element of B. For example, if A = [[a1, a2], [b1, b2], [c1, c2]],
        C = []                             # then B = [ [a1, b1, c1], [a1, b1, c2], [a1, b2, c1], [a1, b2, c2], ...  ]
        for yy in B:
            C += [reset(term, indlist, yy)]
        D = []
        for zz in C:
            if(explicit_spin_check(zz, 1) == 1):
                D += [zz]
        return D       



def explicit_spin_check(term, allow_general_inds):       # make sure in each mat element and operator, there are the same number of alpha and beta indices; 
    result = 1                                           # return 0, if the numbers are different. allow_general_inds determines if we allow some 'general' indices 
    for xx in term.coefficient.matElement:               # first, check mat elements
        if(explicit_spin_check_help(xx.matUpperIndicees, xx.matLowerIndicees, allow_general_inds) == 0):
            result = 0
        if(result == 0):
            break
    if(result == 0):
        return result
    else:                                                # check operator part
        return explicit_spin_check_help(term.uOperators[0].upperIndicees, term.uOperators[0].lowerIndicees, allow_general_inds)



def explicit_spin_check_help(upps, lows, allow_general_inds):    # check if 'upps' and 'lows' have the same nubmer of alpha and beta indices
    if(len(upps) == 0):
        return 1
    else:
        result = 1
        u_alpha = 0
        u_beta  = 0
        u_ab    = 0                                              # number of 'general' (i.e. spin implicit) indices in 'upps'  
        l_alpha = 0
        l_beta  = 0
        l_ab    = 0
        for yy in upps:
            if(yy.other[0] == ""):
                u_ab += 1
            elif(yy.other[0] == 'explicit' and yy.type in ['P', 'E', 'H', 'M']):
                u_alpha += 1
            elif(yy.other[0] == 'explicit' and yy.type in ['p', 'e', 'h', 'm']):
                u_beta += 1
            else:
                print "\n\n                  index.other[0] unexpected: explicit_spin_check_help@real.py\n\n"
                printterm(term)
                print yy.gvIndex
                pdb.set_trace()
        for zz in lows:
            if(zz.other[0] == ""):
                l_ab += 1
            elif(zz.other[0] == 'explicit' and zz.type in ['P', 'E', 'H', 'M']):
                l_alpha += 1
            elif(zz.other[0] == 'explicit' and zz.type in ['p', 'e', 'h', 'm']):
                l_beta += 1
            else:
                print "\n\n                  index.other[0] unexpected: explicit_spin_check_help@real.py\n\n"
                printterm(term)
                print yy.gvIndex
                pdb.set_trace()
        if(not allow_general_inds):
            if(u_ab > 0 or l_ab > 0):
                print "\n general indices not allowed  \n" 
                set_trace()
        if(u_alpha + u_ab >= l_alpha - l_ab and l_alpha + l_ab >= u_alpha - u_ab):
            pass                            # make sure that, in the presence of some 'general' indices, 'upps' and 'lows' won't definitely have different numbers of alpha indices 
        else:
            result = 0
        return result




correspond_spatial_tensor = {'v': 'G',   # when convert spin tensors to tensors with spatial orbitals (antisymmetry is then lost), this dictionary gives the name correspondence
                             'r': 'R'}   # the first one is two-electron repulsion integral; the second one is f12 integral





def convert_spin_to_spatial_tensor_s(terms, spin_tensor_name):              
    result = []
    for xx in terms:
        result +=  convert_spin_to_spatial_tensor(xx, spin_tensor_name)
    return result


def convert_spin_to_spatial_tensor(term, spin_tensor_name):
    counter = 0
    for xx in term.coefficient.matElement:
        if(xx.name == spin_tensor_name):
            counter += 1
    if(counter == 0):
        return [term]
    else:
        result = [copy.deepcopy(term)]                                         
        for i in range(counter):
            result = convert_spin_to_spatial_tensor_1_s(result, spin_tensor_name)   
        return result


def convert_spin_to_spatial_tensor_1_s(terms, spin_tensor_name):
    result = []
    for xx in terms:
        result += convert_spin_to_spatial_tensor_1(xx, spin_tensor_name)
    return result


def convert_spin_to_spatial_tensor_1(term, spin_tensor_name):    # assume that all upper and lower indices are ordered in the same order: either alpha-beta or beta-alpha
    flag = 0
    others = []
    flag2 = 1
    for xx in term.coefficient.matElement:
        if(xx.name == spin_tensor_name and flag2 == 1):           
            flag += 1
            tensor = xx
            flag2 = 0
        else:
            others += [xx]
    if(flag == 0):
        printterm(term)
        set_trace()
        return []
    elif(len(tensor.matUpperIndicees) != 2 or spin_tensor_name not in correspond_spatial_tensor.keys()):
        print "\n not implemented for this tensor\n"
        set_trace()
    else:
        Us = tensor.matUpperIndicees
        Ls = tensor.matLowerIndicees
        upper_letter_range = range(ord('A'), ord('Z') + 1) 
        lower_letter_range = range(ord('a'), ord('z') + 1) 
        if((ord(Us[0].type) in upper_letter_range and ord(Us[1].type) in upper_letter_range) or (ord(Us[0].type) in lower_letter_range and ord(Us[1].type) in lower_letter_range)):
            return [Term(Coefficient(term.coefficient.const, others + [MatElement(correspond_spatial_tensor[spin_tensor_name], Us, Ls)]), term.uOperators),
                    Term(Coefficient(-1.0*term.coefficient.const, others + [MatElement(correspond_spatial_tensor[spin_tensor_name], [Us[1], Us[0]], [Ls[1], Ls[0]])]), term.uOperators)]
        else:
            return [Term(Coefficient(term.coefficient.const, others + [MatElement(correspond_spatial_tensor[spin_tensor_name], Us, Ls)]), term.uOperators)]


        




                  
#--------------------------------------------------------------------------------------------------------------------------------------
                  
#                                                 SRCC

#--------------------------------------------------------------------------------------------------------------------------------------




def srcc_resid_contraction(terms, Ts, order):
    """ Contraction of 'terms' and ex operators 'Ts'; the resulting terms may be further contracted with other ex operators.

        The goals is to obtain contributions for single reference residual equation; since the ops are pure excitation,
        we can discard some terms in the middle. order = 1 for singles and 2 for doubles ... """
    res = []
    binarytuples = [(xx, yy) for xx in terms for yy in Ts]
    #print "total # of terms: to be processed: " + str(len(terms) * len(Ts)) + '\n\n'
    #counter = 0
    for xx, yy in binarytuples:
        #counter += 1
        #if counter % 10 == 0: print "the " + str(counter) + "-th term\n"
        xxrank, yyrank = xx.give_ex_deex_rank(), yy.give_ex_deex_rank()
        if xxrank[0] + yyrank[0] - xxrank[1] > 2 * order: 
            continue  # with or without further contraction, this won't contribute
        if 0:
            print "contribute@srcc_resid_contraction"
            printterms([xx, yy])
        a = contraction_term_term(xx, yy)
        res.extend([zz for zz in a if zz.give_ex_deex_rank()[0] <= 2 * order])
    #print "done\n\n"
    #set_trace()
    return res
    























                  
#--------------------------------------------------------------------------------------------------------------------------------------
                  
#                                                 MR-F12

#--------------------------------------------------------------------------------------------------------------------------------------
        











class MRF12Separator:
    pass








# this section is for the universal perturbation [2]R12 theory of Torhyden and Valeev


def convert_some_general_indices_to_hole(termm):
    """
        convert those 'general indices' coming from resolution of identity used in F12 theory to hole indices, which appear in density or cumulants;
        this function is initially written for the trial Valeev's MR-R12 theory; so far we only have three classes of indices, hole/particle/general,
        which are used to emulate the OBS/external/CBS 
    """
    set_trace()
    indices_to_be_replaced = []
    new_inds = [] 
    for xx in termm.coefficient.matElement:
        if(xx.name in dm_or_cumu_names):
            for yy in xx.matUpperIndicees + xx.matLowerIndicees:
                if(yy.type == 'g'):
                    indices_to_be_replaced += [yy]
    clock = 1
    for xx in indices_to_be_replaced:
        yy = copy.deepcopy(xx)
        yy.type = 'h'
        yy.num = str(1000 + clock)
        new_inds += [yy]
        clock += 1
    termm = reset(termm, indices_to_be_replaced, new_inds)
    termm = dress([termm])[0]
    return termm


            
def convert_some_general_indices_to_hole_s(terms):
    result = []
    for xx in terms: result += [convert_some_general_indices_to_hole(xx)]
    return result
    
    

def convert_some_general_indices_to_particle(termm):
    """
        in the current version of MR-F12 theory of Valeev, at least one of the two creation operators in F12-Gamma-operator (refer to Valeev's first MR-F12 paper)
        must be an 'external' index, which is emulated as a 'particle' index here; for the two indices corresponding to those in the creation operators, if one
        of them has been converted to a 'hole'/OBS index, the other one must be an 'external'/particle index; this function should be used AFTER 
        convert_some_general_indices_to_hole_s 
    """
    print "\n I am a little confused about the implementation of this function; check before use\n"
    set_trace()
    inds_to_be_replaced = []
    new_inds = []
    valid = 1
    for xx in termm.coefficient.matElement:
        if(xx.name == 'r' and valid == 1):
            upps = xx.matUpperIndicees
            lows = xx.matLowerIndicees
            if(upps[0].type == 'h' and upps[1].type == 'h' and lows[0].type == 'h' and lows[1].type == 'h'):
                valid = 0              # this term should be discarded
            elif(upps[0].type == 'h' and upps[1].type == 'h'):
                if(lows[0].type == 'h' and lows[1].type == 'g'): inds_to_be_replaced += [lows[1]]
                elif(lows[1].type == 'h' and lows[0].type == 'g'): inds_to_be_replaced += [lows[0]]
                elif(lows[1].type == 'g' and lows[0].type == 'g'): pass
                else: set_trace()
            elif(lows[0].type == 'h' and lows[1].type == 'h'):
                if(upps[0].type == 'h' and upps[1].type == 'g'): inds_to_be_replaced += [upps[1]]
                elif(upps[1].type == 'h' and upps[0].type == 'g'): inds_to_be_replaced += [upps[0]]
                elif(upps[1].type == 'g' and upps[0].type == 'g'): pass
                else: set_trace()
            else: set_trace()
    if(valid == 0): return []
    else:
        clock = 1
        for xx in indices_to_be_replaced:
            yy = copy.deepcopy(xx)
            yy.type = 'p'
            yy.num = str(1000 + clock)
            new_inds += [yy]
            clock += 1
        termm = reset(termm, inds_to_be_replaced, new_inds)
        termm = dress([termm])[0]
    return termm

        

def convert_some_general_indices_to_particle_s(terms):
    result = []
    for xx in terms: result.extend(convert_some_general_indices_to_particle(xx))
    return result
    


def mr_f12_valeev_discard_s(terms):
    """
        implement the criterion by which some terms are discarded:
        discard terms containing quadratic (or more) 2-body cumulants;
        terms in which two R matrix elements or a R and a g/v are connected via a 2-body cumulants; so far we make the
        following assumption which holds in the intial Valeev theory: every term either contains only one R and only one G or contains only two R without g/v 
    """
    return mr_f12_valeev_discard_s_step_2(mr_f12_valeev_discard_quadra_cumu(terms))




def mr_f12_valeev_discard_quadra_cumu(terms):
    """ 1-body cumulant is called 'lambda1'/'Lambda1'; the maximum rank of kept cumulant is usually 2; so the next line discards quardratic 2-body cumulants, etc """
    if(system_switch['use spin orbital'] == 1):
        terms = select_terms_rank_n(terms, [], ['lambda'], 0) + select_terms_rank_n(terms, ['lambda'], [], 1)              
    else:
        terms = select_terms_rank_n(terms, [], ['Lambda'], 0) + select_terms_rank_n(terms, ['Lambda'], [], 1)              
    #result = select_terms_rank_n(terms, ['f'], [], 1)                          # Fock term does not enter the next selection
    #terms = select_terms_rank_n(terms, ['v'], [], 1)
    return terms


def mr_f12_valeev_discard_s_step_2(terms):
    result = []
    for xx in terms:
        if(mr_f12_valeev_discard(xx) == 0): result.append(xx) 
    return result



def mr_f12_valeev_discard(termm):
    """ this function discards the term in which r/g are DIRECTLY connected via a two-body cumulant """
    if(system_switch['use spin orbital'] == 0): LambdaName = 'Lambda'
    else: LambdaName = 'lambda'
    discard = 0
    count = 0
    L = []
    R = []
    CU = []
    L_indicator = 1
    R_indicator = 1
    for xx in termm.coefficient.matElement:
        if(xx.name in ['r', 'v'] and L_indicator):
            L += [xx]
            L_indicator = 0
        elif(xx.name in ['r'] and R_indicator):
            R += [xx]
            R_indicator = 0
        elif(xx.name in ['r', 'v']):
            print " \n error: more than two R/V \n"
            set_trace()
        elif(xx.name == LambdaName and len(xx.matUpperIndicees) == 2):
            CU += [xx]
    if(len(L) == len(R) == 1):
        pass
    else:
        print "\n something unexpected happens; not necessarily wrong \n"
        printterm(termm)
        set_trace()
    if(len(CU) == 0):
        pass
    elif(len(CU) > 1):
        printterm(termm)
        set_trace() 
    else:
        Linds = L[0].matUpperIndicees + L[0].matLowerIndicees
        Rinds = R[0].matUpperIndicees + R[0].matLowerIndicees  
        CUinds = CU[0].matUpperIndicees + CU[0].matLowerIndicees
        if(overlap(Linds, CUinds) and overlap(Rinds, CUinds)):
            discard = 1    
    return discard  







def pre_screen_2_b_cu(ta, ob, oc):
    """
        This function is used in the next function. 
        The goal is: if the 'lambda' in ta comes from contractions which are not purely LU-RL or LL-RU, then return 1; otherwise 0;
        The current function only implements a prescreening: it tests whether the cumulant is connected to two and only two groups of indices among
        the four groups from 'ob' and 'oc'. 'ob' and 'oc' are operators. 
    """
    if(system_switch['use spin orbital'] == 1): CumuName = 'lambda'
    else: CumuName = 'Lambda'
    result = 0
    lams = []
    for xx in ta.coefficient.matElement:
        if(xx.name == CumuName and len(xx.matUpperIndicees) == 2): lams += [xx]
    if(len(lams) == 0): return result
    elif(len(lams) != 1):
        print " unexpected \n"
        set_trace()
    else:
        inds1 = ob.upperIndicees
        inds2 = ob.lowerIndicees
        inds3 = oc.upperIndicees
        inds4 = oc.lowerIndicees
        overlapp = 0
        yyinds = lams[0].matUpperIndicees + lams[0].matLowerIndicees
        for zz in [inds1, inds2, inds3, inds4]:
            if(overlap(zz, yyinds)): overlapp += 1
        if(overlapp != 2): result = 1
        return result
          
          
               



                                                                                                        
def combine_lambda1_and_inverse_to_delta_s(terms):
    """ 
        the matelement 'x1' will be used to simulate the inverse of the metric matrix which is the same as the 1-particle cumulant/DM; this function
        combines the contracted 'x1' and 'lambda1' (i.e. they have one common index) to a delta;
        the current version works both for spin-free and spin-orbital one. 
    """
    result = []
    for xx in terms: result += [combine_lambda1_and_inverse_to_delta(xx)]
    return result



def combine_lambda1_and_inverse_to_delta(term):
    """
        here we make the assumption that up to 2 'x1' mats exist, which is true for the current theory 
    """ 
    valid = 0
    if system_switch['use spin orbital']:
        MetricMatName = 'lambda1'
    else:
        MetricMatName = 'Lambda1'
    lam = []
    lam_inv = []
    others = []
    for xx in term.coefficient.matElement:
        if(xx.name == 'x1'):
            lam_inv += [xx]
        elif(xx.name == MetricMatName):
            lam += [xx]
        else:
            others += [xx]
    if(len(lam_inv) == 0):
        return term
    elif(len(lam_inv) > 2):
        print " two many 'x1'; either error or the function needs to extended "
        set_trace()
    elif(len(lam_inv) == 1):
        find = 0
        lam1 = []
        lam2 = []
        for yy in lam:
            if not find:
                if yy.matUpperIndicees[0].equal(lam_inv[0].matLowerIndicees[0]):
                    find = 1
                    lam1.append(yy) 
                    continue
                elif yy.matLowerIndicees[0].equal(lam_inv[0].matUpperIndicees[0]):
                    find = 2
                    lam1.append(yy) 
                    continue
                else:
                    lam2.append(yy) 
            else:
                lam2.append(yy)
        if find:
            dicc = {}
            dicc[1] = [MatElement('delta', lam_inv[0].matUpperIndicees, lam1[0].matLowerIndicees)]
            dicc[2] = [MatElement('delta', lam1[0].matUpperIndicees, lam_inv[0].matLowerIndicees)]
            return Term(Coefficient(term.coefficient.const, others + lam2 + dicc[find]), term.uOperators) 
        else:
            print " the mat corresponding to the inverse not found "
            printterm(term)
            #pdb.set_trace() 
            return term
    else:
        lama = []
        lamb = []
        lam2 = []
        find1 = 0
        find2 = 0
        for yy in lam: 
            yymatch = 0
            if find1 == 0:
                if yy.matUpperIndicees[0].equal(lam_inv[0].matLowerIndicees[0]):              
                    find1 = 1
                    lama.append(yy)
                    yymatch = 1
                    continue
                elif yy.matLowerIndicees[0].equal(lam_inv[0].matUpperIndicees[0]):            
                    find1 = 2
                    lama.append(yy) 
                    yymatch = 1
                    continue
            if find2 == 0:
                if yy.matUpperIndicees[0].equal(lam_inv[1].matLowerIndicees[0]):              
                    find2 = 1
                    lamb.append(yy) 
                    yymatch = 1
                    continue
                elif yy.matLowerIndicees[0].equal(lam_inv[1].matUpperIndicees[0]):            
                    find2 = 2
                    lamb.append(yy) 
                    yymatch = 1
                    continue
            lam2.append(yy)                                                                            
        if find1 > 0 and find2 > 0:
            dicc = {}
            dicc[1] = [MatElement('delta', lam_inv[0].matUpperIndicees, lama[0].matLowerIndicees)]
            dicc[2] = [MatElement('delta', lama[0].matUpperIndicees, lam_inv[0].matLowerIndicees)]
            dicc[3] = [MatElement('delta', lam_inv[1].matUpperIndicees, lamb[0].matLowerIndicees)]
            dicc[4] = [MatElement('delta', lamb[0].matUpperIndicees, lam_inv[1].matLowerIndicees)]
            return Term(Coefficient(term.coefficient.const, others + lam2 + dicc[find1] + dicc[2 + find2]), term.uOperators)    
        elif find1 > 0: 
            dicc = {}
            dicc[1] = [MatElement('delta', lam_inv[0].matUpperIndicees, lama[0].matLowerIndicees)]
            dicc[2] = [MatElement('delta', lama[0].matUpperIndicees, lam_inv[0].matLowerIndicees)]
            return Term(Coefficient(term.coefficient.const, others + lam2 + [lam_inv[1]] +  dicc[find1]), term.uOperators)    
        elif find2 > 0:
            dicc = {}
            dicc[1] = [MatElement('delta', lam_inv[1].matUpperIndicees, lamb[0].matLowerIndicees)]
            dicc[2] = [MatElement('delta', lamb[0].matUpperIndicees, lam_inv[1].matLowerIndicees)]
            return Term(Coefficient(term.coefficient.const, others + lam2 + [lam_inv[0]] + dicc[find2]), term.uOperators)    
        else:
            print " the mat corresponding to the inverse not found "                                    
            printterm(term)                                                                             
            #pdb.set_trace()                                                                             
            return term
        


def group_mats_by_indices_s(terms, extinds, deal):
    """ 
        for every term, reorder its mats: firstly, those mats which has indices in 'extinds' or which connects those mats; then the other mats.
        This function is not general, and it only works for the current theory 
    """
    result = []
    for xx in terms: result += group_mats_by_indices(xx, extinds, deal)
    return result



def group_mats_by_indices(term, extinds, deal):
    inds_basis = copy.deepcopy(extinds)
    mats_r = []  # mats containing 'extinds'
    mats_o = []  # others
    matsc = []   # mats in 'others' which are connected to mats_r 
    matso = []   # other mats in 'others'
    for xx in term.coefficient.matElement:
        if(xx.name == 'r'):            # the f12 terms contains extinds
            mats_r += [xx]
            inds_basis += xx.matUpperIndicees + xx.matLowerIndicees
        else: mats_o += [xx]
    for yy in mats_o:
        if(overlap(yy.matUpperIndicees + yy.matLowerIndicees, inds_basis)):
            matsc += [yy]              # quantities connecting two f12 terms 
            #printMatElements(matsc)
        else: matso += [yy]
    res = [dress_a(Term(Coefficient(term.coefficient.const, mats_r + matsc + matso), term.uOperators))] # just reorder the mats in 'term'
    if(deal == 'keep'): return res
    else:
        check_ex = 1
        valid = 1
        if(check_ex):
            indss = []
            for zz in mats_r + matsc: indss += zz.matUpperIndicees + zz.matLowerIndicees
            exindss = find_ex_lists(indss, [])   # get the external indices of the mats mats_r + matsc
            if(len(exindss) != len(extinds)): valid = 0
        if(deal == 'del'):
            if(valid): return res
            else: return []
        elif(deal == 'keep_del'):
            if(valid): return []
            else: return res
            





def extract_unique_terms_from_permu_and_conju(terms, pqrs, vwtu):
    
    """ our expression can be written as (X(p,q,r,s;v,w,t,u) + conju(X(v,w,t,u;p,q,r,s)). 
        this function tries to find X so as to reduced the # of terms; for any given term, either the permuted form is equivalant to itself or it is another term in 'terms'
        this function should be used AFTER merge_all.
        The assumption is made that all terms are numbers;i.e., that are no operators 
    """   
    result = []
    result2 = []
    B = copy.deepcopy(terms)
    nn = len(terms)
    count = 0
    while (len(B) > 0):
        if(count == nn):                                                   # the maximum number of cycles is len(B)
            set_trace()
            return
        Y = varyco(-1.0, [invert_upp_low(reset(B[0], pqrs + vwtu, vwtu + pqrs))])
        if(len(merge_all([B[0]] + Y)) == 0):                                  
            result += varyco(0.5, B[0:1])
            del B[0]
        else:
            findd = 0
            lenn = len(B)
            for xx in range(lenn)[1:]:
                if(len(merge_all([ B[xx] ] + Y)) == 0):
                    findd = 1
                    location = xx
                    break
            if(findd == 0):
                printterms(Y)
                #printterms(terms) 
                set_trace()                                                # not explicitly antisymmetrized?
                result += [B[0]]
                del B[0]
            else:
                result += [B[0]]
                del B[location]
                del B[0]
        count += 1
    return result



def permu_and_conju_s(terms, pqrs, vwtu):
    """ the reverse function of the above """

    result = []
    for xx in terms: result += [canonical_one_by_permu(invert_upp_low(reset(xx, pqrs, vwtu)))] 
    return result

    

def replace_sf_3_dm_with_cumu_s(terms):              
    """ replace spin-free 3-RDM by 1/2-cumulant (spin-free), as in the MolPhys paper by Kutzelnigg, Eq. (75)  """

    result = []
    for xx in terms: result += replace_sf_3_dm_with_cumu(xx)
    return result



def replace_sf_3_dm_with_cumu(term):
    counter = 0
    for xx in term.coefficient.matElement:
        if(xx.name == 'd3'):
            counter += 1
    if(counter == 0): return [term]
    else:
        def replace_sf_3_dm_with_cumu_1_s(terms):
            result = []
            for xx in terms: result.extend(xx.replace_sf_3_dm_with_cumu_1())
            return result
        result = [copy.deepcopy(term)]                     
        for i in range(counter): result = replace_sf_3_dm_with_cumu_1_s(result)        
        return result












def find_pair_pp_indices(termm):
    """ in termm, if a geminal element contains two upper or lower CABS indices, return them; assume there are at most always 2 CABS indices """
    result = []
    for xx in termm.coefficient.matElement:
        if(xx.name == 'r'):
            linds = xx.matLowerIndicees
            uinds = xx.matUpperIndicees
            if(linds[0].type == linds[1].type == 'p'):
                result += linds
            elif(uinds[0].type == uinds[1].type == 'p'):
                result += uinds
    return uniqueinds_from_indlist(result)



""" for [2]r12, there are at most 3 CABS indices;
    to remove mats of 2 CABS indices, I will first replace 'p' by 'g' and 'h'; then classify terms into 3 groups: 1-g, 2-g, 3-g;
    1 keep 3-g unchanged; 
    2 break 1-g; 
    3 for 2-g,
      if two 'r' (metric matrix) are contracted by the two 'g' indices, keep unchanged;
      if one 'r' (r1) has 2-g and then contracted to 'f' and the other 'r' (r2), break the 'g' in r2; 
         then, for the broken terms, there is one 'g' left 
         for the last 'g' index in 'r1', (1) if its neighbor is 'h', then break 
                                         (2) if its neighbor is 'p', then keep 'g'
      if each 'r' has only one 'g', break it 
    The assumption is made that at most 3 CABS indices; all CABS indices appear in 'r' or 'f'
 """


def valeev_f12_resolve_2cabs_terms(terms):
    res = []
    for xx in terms: res.extend(valeev_f12_resolve_2cabs_term(xx))
    return res


def valeev_f12_resolve_2cabs_term(termm):
    res = []
    uniqueinds = give_indlist(termm)
    if 0:
        printterm(termm) 
        showindlist(uniqueinds)
    ps = pick_type_ind_from_inds(uniqueinds, 'p') # replace CABS indices by CABS and OBS
    replace_ps = expand_specific_ps(termm, ps)
    if 0:
        set_trace()
        printterms(replace_ps)
    for yy in replace_ps:
        if 0: set_trace()
        gt = valeev_f12_resolve_2cabs_deal_with_g_term(yy)
        if 0:
            set_trace()
            printterm(yy) 
            print "\n\n"
            printterms(gt) 
        res.extend(gt) 
    return res


def valeev_f12_resolve_2cabs_deal_with_g_term(termm):
    """ try to determine whether or not and how to expand 'g' indices """
    #set_trace()
    uniqueinds = give_indlist(termm)
    numtype = counttype(uniqueinds)
    if numtype[4] > 3:
        print "\n should not happen \n"; set_trace()
    if numtype[4] in [0, 3]:
        return [termm]
    if numtype[4] in [1]:
        return break_g(termm)
    else:
        rs = []
        for xx in termm.coefficient.matElement:
            if xx.name == 'r': rs.append(xx)
            elif xx.name == 'f': f = xx
        if len(rs) == 1: return [termm]    # this is the linear term < V \Omega_1 >; this term has been dealt with properly so far by default, since terms of 2 'g' indices 
        r1 = rs[0]                         # should just stay unchanged for this term
        r2 = rs[1]
        r1inds = r1.matUpperIndicees + r1.matLowerIndicees
        r2inds = r2.matUpperIndicees + r2.matLowerIndicees
        numtype1 = counttype(r1inds)
        numtype2 = counttype(r2inds)
        if numtype1[4] == numtype2[4] == 2: return [termm] # 2 CABS indices are contracted betwen two geminal
        elif numtype1[4] == numtype2[4] == 1: return break_g(termm)
        else:
            res = []
            if numtype1[4] == 2 and numtype2[4] == 1:
                oneg = pick_type_ind_from_inds(r2inds, 'g')
                temps = expand_specific_g(termm, oneg)
            elif numtype1[4] == 1 and numtype2[4] == 2:
                oneg = pick_type_ind_from_inds(r1inds, 'g')
                temps = expand_specific_g(termm, oneg)
            else: set_trace()
            for xy in temps:
                res += valeev_f12_special_treat_1_g(xy)
            return res



def valeev_f12_special_treat_1_g(termm):
    """ this is a special treatment at the end of the whole procedure"""
    if 0:
        set_trace()
        printterm(termm) 
    rs = []
    for xx in termm.coefficient.matElement:
        if xx.name == 'r': rs.append(xx)
    r1 = rs[0]
    r2 = rs[1]
    r1inds = rs[0].matUpperIndicees + rs[0].matLowerIndicees
    r2inds = rs[1].matUpperIndicees + rs[1].matLowerIndicees
    numg1 = counttype(r1inds)[4]
    numg2 = counttype(r2inds)[4]
    if numg1 == 1: ma = r1inds
    else: ma = r2inds
    g1 = pick_type_ind_from_inds(ma, 'g')[0] 
    pos = g1.gvIndexNo(ma)
    neighbor = {0:1, 1:0, 2:3, 3:2}[pos]
    if ma[neighbor].type == 'h': return break_g(termm)
    elif ma[neighbor].type == 'p': return [termm]
    else: set_trace()
    
           

    
                
            
             

#--------------------------------------------------------------------------------------------------------------------------------------

#                                                 Factorization

#--------------------------------------------------------------------------------------------------------------------------------------
# Factorization related functions






class FactorizationSeparator:
    pass





def group_factorize_terms_based_on_mats_list(terms, mats_list):
    """ sequential factorization """
    fac = [[terms]]
    for xx in mats_list:
        res = separate_terms_wrt_factorization(fac[-1][0], xx)
        fac = fac[:-1] + res
    return fac

    

    
def separate_terms_wrt_factorization(terms, mats):
    """seperate terms into 2 groups: those containing mats and those not. if terms = [0.5 AB, AC, DE] and mats = [A], return [[[A], [0.5B, C]], [DE]]"""
    TS1 = []
    TS2 = []
    for xx in terms:
        inf = xx.hasmats(mats)
        if inf[0]: TS1.append(Term(Coefficient(xx.coefficient.const, inf[1]), xx.uOperators, xx.type, xx.other))
        else: TS2.append(xx)
    if 0:
        printterms(terms)
        set_trace()
    return [[mats, TS1], [TS2]] 


def group_terms_based_on_mats_list(terms, mats_list):
    """ basically, this function groups terms, firstly based on whether a terms contains mats or not; then for the groups terms, first classify each group according to ... """
    res = [terms]
    result = []
    for xx in mats_list: res =  group_termslist_based_on_mats(res, xx)
    for yy in res: result.extend(yy)
    return result



def group_termslist_based_on_mats(termslist, mats):
    res = []
    for xx in termslist: res.extend(group_terms_based_on_mats(xx, mats))
    return res


def group_terms_based_on_mats(terms, mats):
    """group terms into 2 groups"""
    TS1 = []
    TS2 = []
    for xx in terms:
        inf = xx.hasmats(mats)
        if inf[0]: TS1.append(xx)
        else: TS2.append(xx)
    return [TS1, TS2] 



def group_terms_based_on_matname_list(terms, matname_list):
    """ basically, this function groups terms, firstly based on whether a terms contains a mat named 'xx'; then for the groups terms, first classify each group according to ... """
    res = [terms]
    result = []
    for xx in matname_list: res = group_termslist_based_on_matname(res, xx)
    for yy in res: result.extend(yy)
    return result



def group_termslist_based_on_matname(termslist, matname):
    res = []
    for xx in termslist: res.extend(group_terms_based_on_matname(xx, matname))
    return res


def group_terms_based_on_matname(terms, matname):
    """ return a list of two groups, depending on the term has a mat named 'matname' or not; no effect on an individual term """
    t1 = []
    t2 = []
    for xx in terms:
        if xx.hasmat_name(matname)[0]: t1.append(xx)
        else: t2.append(xx)
    return [t1, t2]



        




#--------------------------------------------------------------------------------------------------------------------------------------

#                                                 LATEX

#--------------------------------------------------------------------------------------------------------------------------------------







class LatexSeparator:
    pass









def termslist_texform(sf_or_so, termslistt, filename, EqnName = '', comment = "", WhichMap = 'sf:valeev:mrf12:new', terms_one_line = 3, readable_form = 1, print_total_num = 0, beginarray =1, endarray = 1):          
    """  
         termslist is a list of lists, [ [mats, other contribution], [], ...]; this function is for (partially) factorized output;
         'mat' is the common factor; in principle, more general factorization should be allowed, but it may suffice for now.
         sf_or_so == 'sf' or 'so';
         filename is the name of the latex output file;
         EqnName is the name of the left value in the eqnarray, e.g. A = ....;
         comment is some sentence before the equations; 
         terms_one_line is the number of terms on each line;
         if 'readable_form' == 1, we would replace p1/p2 by a/b, etc; 'readable_form' does not apply to indices in operators!!
    """ 
    if(1):                                                                     
        termslist= termslistt
        res = "\documentclass[a4paper, 12pt]{article}\n\\usepackage{latexsym, amsfonts}\n \\begin{document}\n"
        res += comment + "\n"
        if(beginarray):
            res += "\\begin{eqnarray}\n"
        totalnum = len(termslist)
        counter = 0
        for xx in termslist:
            if(counter == 0):
                res += "%s & = & " % EqnName
            elif len(xx) == 0 or (len(xx) == 1 and not xx[0]) or (len(xx) == 2 and not xx[1]): continue
            else: res += " &  &  " 
            if(len(xx) == 2): # with the common factor
                num_commonmats = len(xx[0])
                res += texform_mats(xx[0], readable_form, WhichMap)
                terms2 = xx[1] 
            elif(len(xx) == 1): # without the common factor
                num_commonmats = 0
                terms2 = xx[0] 
            else:
                print len(xx)
                set_trace()
            res += "\\big( "
            if(1):
                for i in range(len(terms2)):
                    termconst = terms2[i].coefficient.const 
                    res += texform_const_coeff(termconst)
                    termmats = terms2[i].coefficient.matElement
                    for k in range(len(termmats)):
                        matk = termmats[k]
                        res += texform_one_mat(matk, readable_form, WhichMap)
                    OP = terms2[i].uOperators[0]
                    res += texform_one_operator(OP, readable_form, WhichMap)
                    #if num_commonmats == 2 and i == 2: set_trace()
                    if(counter == totalnum-1):
                        res += termslist_texform_get_end(num_commonmats, i, terms2, terms_one_line, 1)
                    else:
                        res += termslist_texform_get_end(num_commonmats, i, terms2, terms_one_line, 0)
            counter += 1
        if(endarray):
            res += "\end{eqnarray}\n"
        if(print_total_num):
            res += "No. of terms = %s \\\\\n" % len(terms2)
        res += "\end{document}"
        file = open(filename, "w")
        file.write(res)
        file.close()







def texform(sf_or_so, terms, filename, EqnName = '', comment = "", WhichMap = 'sf:valeev:mrf12:new', terms_one_line = 3, readable_form = 1, print_total_num = 0):          
    """  sf_or_so == 'sf' or 'so';
         filename is the name of the latex output file;
         EqnName is the name of the left value in the eqnarray, e.g. A = ....;
         comment is some sentence before the equations; 
         terms_one_line is the number of terms on each line;
         if 'readable_form' == 1, we would replace p1/p2 by a/b, etc; 'readable_form' does not apply to indices in operators!!
         
    """ 
    if(1):                                                                     
        res = "\documentclass[a4paper, 12pt]{article}\n\\usepackage{latexsym, amsfonts}\n \\begin{document}\n"
        res += texform_mainbody(sf_or_so, EqnName, comment, terms, terms_one_line, 1, 1, readable_form, WhichMap, print_total_num)     
        res += "\end{document}"
        file = open(filename, "w")
        file.write(res)
        file.close()





def texform_mainbody(sf_or_so, EqnName, comment, terms2, terms_one_line, beginarray, endarray, readable_form, WhichMap, print_total_num):
    if(1):
        res = comment + "\n"
        if(beginarray):
            res += "\\begin{eqnarray*}\n"
        res += "%s & = & " % EqnName
        for i in range(len(terms2)):
            termconst = terms2[i].coefficient.const 
            res += texform_const_coeff(termconst)
            termmats = terms2[i].coefficient.matElement
            for k in range(len(termmats)):
                matk = termmats[k]
                res += texform_one_mat(matk, readable_form, WhichMap)
            OP = terms2[i].uOperators[0]
            res += texform_one_operator(OP, readable_form, WhichMap)
            res += texform_get_end(i, terms2, terms_one_line)
        if(endarray):
            res += "\end{eqnarray*}\n"
        if(print_total_num):
            res += "No. of terms = %s \\\\\n" % len(terms2)
        return res




def texform_const_coeff(termconst):
    if(1):
        if(1):
            res = ""
            if(termconst == 1.0):
                res += "\ "
            elif(termconst == -1.0):
                res += "-"
            elif(abs(termconst - 0.5) < 1e-5):
                res += "\\frac{1}{2} "
            elif(abs(termconst + 0.5) < 1e-5):
                res += "-\\frac{1}{2} "
            elif(abs(termconst - 1.0/3.0) < 1e-5):
                res += "\\frac{1}{3} "
            elif(abs(termconst + 1.0/3.0) < 1e-5):
                res += "-\\frac{1}{3} "
            elif(abs(termconst - 0.25) < 1e-5):
                res += "\\frac{1}{4} "
            elif(abs(termconst + 0.25) < 1e-5):
                res += "-\\frac{1}{4} "
            elif(abs(termconst - 0.2) < 1e-5):
                res += "\\frac{1}{5} "
            elif(abs(termconst + 0.2) < 1e-5):
                res += "-\\frac{1}{5} "
            elif(abs(termconst - 1/8.0) < 1e-5):
                res += "\\frac{1}{8} "
            elif(abs(termconst + 1/8.0) < 1e-5):
                res += "-\\frac{1}{8} "
            elif(termconst > 0):
                res += str(termconst)+ " "
            else:
                res += "\! " + str(termconst)+ " "
            return res


def texform_mats(mats, readable_form, WhichMap):
    ss = ""
    for matk in mats: ss += texform_one_mat(matk, readable_form, WhichMap)
    return ss

def texform_one_mat(matk, readable_form, WhichMap):
    if(1):
        if(1):
            if(1):
                res = dic_matname_va[matk.name] + "^{"
                for m in range(len(matk.matUpperIndicees)):
                    if(readable_form):
                        res += convert_to_readable_form(WhichMap, matk.name, matk.matUpperIndicees[m].type, matk.matUpperIndicees[m].num)
                    else:
                        res += matk.matUpperIndicees[m].type + "_{" + matk.matUpperIndicees[m].num + "} "
                res = res + "}_{"
                for m in range(len(matk.matLowerIndicees)):
                    if(readable_form):
                        res += convert_to_readable_form(WhichMap, matk.name, matk.matLowerIndicees[m].type, matk.matLowerIndicees[m].num)
                    else:
                        res += matk.matLowerIndicees[m].type + "_{" + matk.matLowerIndicees[m].num + "} "
                res += "}  "
                return res






def texform_one_operator(OP, readable_form, WhichMap):
    if(1):
        if(1):
            res = ""
            if(len(OP.upperIndicees) == 0):
                pass
            else:
                if not system_switch['use spin orbital']:
                    res += "\hat{e}^{"
                else:
                    res += "a^{"
                for j in range(len(OP.upperIndicees)):
                    if readable_form:
                        res += convert_to_readable_form(WhichMap, "", OP.upperIndicees[j].type, OP.upperIndicees[j].num)
                    else:
                        res += OP.upperIndicees[j].type + "_{" + OP.upperIndicees[j].num + "} "
                res += "}_{"
                for j in range(len(OP.lowerIndicees)):
                    if readable_form:
                        res += convert_to_readable_form(WhichMap, "", OP.lowerIndicees[j].type, OP.lowerIndicees[j].num)
                    else:
                        res += OP.lowerIndicees[j].type + "_{" + OP.lowerIndicees[j].num + "} "
                res +=  "}"
            return res    



def texform_get_end(i, terms2, terms_one_line):
    if(1):
        if(1):
            res = ""
            if terms_one_line == 1:
                if i < len(terms2)-1:
                    if terms2[i+1].coefficient.const < 0.0:
                        res += "\\nonumber \\\\\n& &"
                    else:
                        res += "+\\nonumber \\\\\n& & "
                else:
                    res +=  "\\nonumber \n"
            elif i == 0:
                if i < len(terms2)-1:
                    if terms2[i+1].coefficient.const < 0.0:
                        res += " "
                    else:
                        res += "+ "
                else:
                    res +=  '\\nonumber \n'
            elif i % terms_one_line == 0:
                if i < len(terms2)-1:
                    if(terms2[i+1].coefficient.const < 0.0):
                        res += " "
                    else:
                        res += "+ "
                else:
                    res +=  ''
                    res +=  "\\nonumber \n"
            elif i % terms_one_line == terms_one_line - 1:
                if i < len(terms2)-1:
                    if terms2[i+1].coefficient.const < 0.0:
                        res += "\\nonumber \\\\\n& &"
                    else:
                        res += "+\\nonumber \\\\\n& & "
                else:
                    res +=  "\\nonumber \n"
            else:
                if i < len(terms2)-1:
                    if terms2[i+1].coefficient.const < 0.0:
                        res += " "
                    else:
                        res += "+ "
                else:
                    res +=  "\\nonumber   \n"
            return res 






def termslist_texform_get_end(num_commonmats, i, terms2, terms_one_line, endoflist):
    if(1):
        if(1):
            res = ""
            space = ""
            one_mat_space = "\quad "
            for x in range(num_commonmats): space += one_mat_space 
            if(terms_one_line == 1):
                if(i < len(terms2)-1):
                    if(terms2[i+1].coefficient.const < 0.0):
                        res += "\\nonumber \\\\\n& &"
                    else:
                        res += "+\\nonumber \\\\\n& & "
                elif(endoflist):
                    res +=  "\\big) \\nonumber  \n"
                else:
                    res +=  "\\big)+ \\nonumber \\\\ \n "
            elif((i - terms_one_line*(i/terms_one_line)) == terms_one_line - 1):
                if(i < len(terms2)-1):
                    if(terms2[i+1].coefficient.const < 0.0):
                        res += "\\nonumber \\\\ \n& &" + space
                    else:
                        res += "+\\nonumber \\\\\n& & " + space
                elif(endoflist):
                    res +=  "\\big) \\nonumber  \n"
                else:
                    res +=  "\\big)+ \\nonumber \\\\ \n "
            else:
                if(i < len(terms2)-1):
                    if(terms2[i+1].coefficient.const < 0.0):
                        res += " "
                    else:
                        res += "+ "
                elif(endoflist):
                    res +=  "\\big) \\nonumber  \n"
                else:
                    res +=  "\\big)+ \\nonumber\\\\ \n "
            return res 





# why these two are special? forgot ...
spatial_index_tensors = ['G', 'R']

def convert_to_readable_form(WhichMap, matname, ttype, nnum):                 # used in texform function to convert index to readable form
    readable_form_dic = ReadableFormMapping[WhichMap]                         # assume in each term, indices like 'p1' 'P1' will not coexist!! 
    keyy = readable_form_dic.keys()
    if(system_switch['use spin orbital'] == 0):        
        if((ttype + nnum) in keyy):
            return readable_form_dic[ttype.lower() + nnum]     
        else:
            return ttype + "_{" + nnum + "} "
    else:
        if(system_switch['explicit spin'] == 0):
            if((ttype + nnum) in keyy):
                return readable_form_dic[ttype.lower() + nnum]     
            else:
                return ttype + "_{" + nnum + "}"
        else:
            print "\n Forget about the use of the following section \n"
            set_trace()
            if(matname in spatial_index_tensors):
                if((ttype + nnum) in keyy):
                    return readable_form_dic[ttype.lower() + nnum]
                else:
                    return ttype + "_{" + nnum + "}"
            else:
                upper_letter_range = range(ord('A'), ord('Z') + 1) 
                lower_letter_range = range(ord('a'), ord('z') + 1) 
                if(ord(ttype) in upper_letter_range):
                    if((ttype + nnum) in keyy):
                        return '{' + readable_form_dic[ttype.lower() + nnum] + '}_{\\alpha}'
                    else:
                        return ttype + "_{" + nnum + "\\alpha}"
                elif(ord(ttype) in lower_letter_range):
                    if((ttype + nnum) in keyy):
                        return '{' + readable_form_dic[ttype.lower() + nnum] + '}_{\\beta}'
                    else:
                        return ttype + "_{" + nnum + "\\beta}"
                else:
                    print matname + ',' + ttype + ',' + nnum + '\n'
                    set_trace()
                    return

 
                     


            

ReadableFormMapping = { 1: {'p1': 'p',
                            'p2': 'q',
                            'p3': 'r',
                            'p4': 's'},
        'sf:valeev:mrf12': {'h1': 'r',
                            'h2': 's',
                            'h3': 'p',
                            'h4': 'q',
                            'h5': 'v',
                            'h6': 'w',
                            'h7': 't',
                            'h8': 'u',
                            'h9': 'x',
                           'h10': 'y',
                           'h11': 'z',
                            'p1': "\\alpha'",
                            'p2': "\\beta'",
                            'p3': "{\\alpha}_2'",
                            'p4': "{\\beta}_2'",
                            'g1': "\kappa",
                            'g2': "\lambda", 
                            'g3': "\mu"},
                   'srsd': {'h1': 'i',
                            'h2': 'j',
                            'h3': 'k',
                            'h4': 'l',
                            'p1': 'a',
                            'p2': 'b',
                            'p3': 'c',
                            'p4': 'd'},
    'sf:valeev:mrf12:new': {'h1': 'p',
                            'h2': 'q',
                            'h3': 'r',
                            'h4': 's',
                            'h5': 'v',
                            'h6': 'w',
                            'h7': 't',
                            'h8': 'u',
                            'h9': 'x',
                           'h10': 'y',
                           'h11': 'z',
                            'p1': "\\alpha'",
                            'p2': "\\beta'",
                            'p3': "{\\alpha}_2'",
                            'p4': "{\\beta}_2'",
                            'g1': "\kappa",
                            'g2': "\lambda", 
                            'g3': "\mu"},
  'sf:valeev:mrf12:verify': {'h1': 'p',
                            'h2': 'q',
                            'h3': 'r',
                            'h4': 's',
                            'h5': 'v',
                            'h6': 'w',
                            'h7': 't',
                            'h8': 'u',
                            'h9': 'x',
                           'h10': 'y',
                           'h11': 'z',
                            'p1': "\\alpha'",
                            'p2': "\\beta'",
                            'p3': "{\\alpha}_2'",
                            'p4': "{\\beta}_2'",
                            'g1': "\mu",
                            'g2': "\\nu", 
                            'g3': "\kappa",
                            'g4': "\lambda"}
                      }












#------------------------------------------------------------------------------------------------------------------------------------------------

#                                                Ohio  Tensor Contraction Engine Interface

#------------------------------------------------------------------------------------------------------------------------------------------------





dic_p_h_2_o_v = {  'v': {'upper': {'p': 'v',
                                   'h': 'o',
                                   'm': 'm',
                                   'e': 'e'},
                         'lower': {'p': 'v',
                                   'h': 'o',
                                   'm': 'm',
                                   'e': 'e'}},
                   'f': {'upper': {'p': 'v',
                                   'h': 'o',
                                   'm': 'm',
                                   'e': 'e'},
                         'lower': {'p': 'v',
                                   'h': 'o',
                                   'm': 'm',
                                   'e': 'e'}},
                  't1': {'upper': {'p': 'v',
                                   'e': 'e'},
                         'lower': {'h': 'o',
                                   'm': 'm'}},
                  't2': {'upper': {'p': 'v',
                                   'e': 'e'},
                         'lower': {'h': 'o',
                                   'm': 'm'}},
                  't3': {'upper': {'p': 'v',
                                   'e': 'e'},
                         'lower': {'h': 'o',
                                   'm': 'm'}},
                   'd': {'upper': {'h': 'o',
                                   'm': 'm'},
                         'lower': {'h': 'o',
                                   'm': 'm'}},
                  's4': {'upper': {'p': 'v',
                                   'e': 'e',
                                   'h': 'h',
                                   'm': 'm'},
                         'lower': {'h': 'o',
                                   'm': 'm'}}}                                 # I should include 's' mats when necessary   




def PHME_2_phme(typee):
    """ Replace upper case by lower case. """
    return typee.replace('P', 'p').replace('H', 'h').replace('M', 'm').replace('E','e')



def phme_2_PHME(typee):
    """ Replace lower case by upper case. """
    return typee.replace('o', 'O').replace('v', 'V').replace('m', 'M').replace('e', 'E')



def find_index_s_type(terms, ttype):
    """ We find all unique indices from 'terms' of type 'ttype' """

    inds = []
    LA = []
    LB = []
    for xx in terms:
        for yy in xx.coefficient.matElement:
            for zz in yy.matUpperIndicees + yy.matLowerIndicees:
                if(zz.type == ttype): inds.append(zz)
    return  combineindlist([], inds)





def mat_spin_state(mat, declare):                                                
    """ Give the spin-state, 'a' means alpha, 'b' means beta; this function is used in next func. When decare == 'implicit', return ''. """

    map_indtype_spin = {'P': 'a', 'H': 'a', 'E': 'a', 'M': 'a', 'p': 'b', 'h': 'b', 'e': 'b', 'm': 'b'}
    if declare == 'explicit':                                                  
        if len(mat.matUpperIndicees) in [1, 2, 3]:
            spinstate = ""
            for xx in mat.matUpperIndicees: spinstate += map_indtype_spin[xx.type]
            return spinstate
        else:
            print "\n\n                    Three-body quantity? not implemented yet.\n"
            set_trace()
    else:
        return ""

    



def TCE_mat_change_mat(mat, declare): 
    """ Return the mat.name in TCE format; also interchange upper and lower inds. """
                                      
    Mat = MatElement('', mat.matLowerIndicees, mat.matUpperIndicees, mat.connectivity, mat.other)
    Name = mat.name[0] + mat_spin_state(mat, declare) + '_'
    for xx in Mat.matUpperIndicees:
        Name += dic_p_h_2_o_v[mat.name]['upper'][PHME_2_phme(xx.type)]
    for xx in Mat.matLowerIndicees:
        Name += dic_p_h_2_o_v[mat.name]['lower'][PHME_2_phme(xx.type)]
    Mat.name = Name
    return Mat



def TCE_term_change_mat(term, declare):                                    
    """ Change all mats from TCE format in the sense of mat.name and the order of upper and lower indices in each mat. """

    coef = Coefficient(term.coefficient.const, [])            
    for xx in term.coefficient.matElement: coef.matElement += [TCE_mat_change_mat(xx, declare)]
    return Term(coef, term.uOperators, term.type, term.other)



def TCE_terms_change_mat(terms, declare):                    
    result = []                                             
    for xx in terms: result.append(TCE_term_change_mat(xx, declare))
    return result




def TCE_get_var_names(terms):                                            
    """ return the unique mat names after doing TCE_term_change_mat(); this is to get the header ncc(in fa_vo[Va, Oa]) """

    result = []
    for xx in terms:
        for yy in xx.coefficient.matElement:
            if yy.name not in result: result.append(yy.name)                      # collect all unique mat names
    result2 = []
    for zz in result:
        print zz                                                               # when this function is done, all the mat names have the  form v_oovv (spatial) or vbb_oovv(spin)
        result2 += [zz + TCE_get_var_names_help(zz[1: zz.index('_')], zz[zz.index('_') + 1:])]
    return result2


def TCE_get_var_names_help(spin, PH):
    S = '['             
    if(spin == ''):  # implicit spin, return [O, V, etc]
        for xx in PH: S += phme_2_PHME(xx) + ","
        S = S[:-1]                                                             # get rid of the last ","   
    elif(spin in ['a', 'b']):
        S += phme_2_PHME(PH[0])+ spin[0] + "," + phme_2_PHME(PH[1]) + spin[0]
    elif(spin in ['aa', 'bb', 'ab']):                                          # return [Oa, Ob,...]
        S += phme_2_PHME(PH[0])+ spin[0] + "," + phme_2_PHME(PH[1]) + spin[1] + "," + phme_2_PHME(PH[2]) + spin[0] + ","  + phme_2_PHME(PH[3]) + spin[1]
    elif(spin in ['aaa', 'aab', 'abb', 'bbb']):                                          # return [Oa, Ob,...]
        S += phme_2_PHME(PH[0])+ spin[0] + "," + phme_2_PHME(PH[1]) + spin[1] + "," \
             + phme_2_PHME(PH[2]) + spin[2] + ","  + phme_2_PHME(PH[3]) + spin[0]  + "," \
             + phme_2_PHME(PH[4]) + spin[1] + ","  + phme_2_PHME(PH[5]) + spin[2]  
    else:
        print "\n\n                    Not implemented yet.\n\n"
        set_trace()
    S += ']'
    return S



def TCE_output_one_term(term):                                                 # output one term A*B*C
    S = ""
    inds = []
    for xx in term.coefficient.matElement:
        inds += xx.matUpperIndicees + xx.matLowerIndicees
    contracted_inds = find_contracted_inds(inds)
    contracted_inds = reorder_inds_hmpe(contracted_inds)
    if term.coefficient.const > 0:
        S += '+' + str(term.coefficient.const)
    elif term.coefficient.const < 0:
        S += str(term.coefficient.const)
    else:
        print "\n\n                                       Zero coefficient? TCE_output_one_term@real.py\n\n"
        set_trace()
    if len(term.coefficient.matElement) == 0:
        print "\n\n                                       No matelement?    TCE_output_one_term@real.py\n\n"
        set_trace()
    elif len(contracted_inds) == 0:
        for mat in term.coefficient.matElement:
            S += "*" + mat.name + '['
            clock = 0
            for indd in mat.matUpperIndicees + mat.matLowerIndicees:
                if clock == len(mat.matUpperIndicees + mat.matLowerIndicees) -1:
                    S += indd.gvIndex()    
                else:
                    S += indd.gvIndex()+ ','
                clock += 1
            S += ']'
    else:
        S += '*sum['
        clock1 = 0
        for mat in term.coefficient.matElement:
            S += mat.name + '['
            clock2 = 0
            for indd in mat.matUpperIndicees + mat.matLowerIndicees:
                if clock2 == len(mat.matUpperIndicees + mat.matLowerIndicees) -1:
                    S += indd.gvIndex()
                else:
                    S += indd.gvIndex()+ ','
                clock2 += 1
            S += ']'
            if clock1 == len(term.coefficient.matElement)-1:
                S += ','
            else:
                S += '*'
            clock1 += 1
        S += '{'
        clock3 = 0
        for contracted_ind in contracted_inds:
            if(system_switch['ccsd(t)'] == 1):
                if(contracted_ind.gvIndex() not in ['p1', 'p2', 'p3', 'h1', 'h2','h3']):
                    S += contracted_ind.gvIndex()
                    if(clock3 != len(contracted_inds)-1):
                        S += ',' 
                    else:
                        pass
                else:
                    if(clock3 != len(contracted_inds)-1):
                        S += ','
                    else:
                        pass
            else:
                S += contracted_ind.gvIndex()
                if(clock3 != len(contracted_inds)-1):
                    S += ','
                else:
                    pass
            clock3 += 1
        if(S[-1] == ','):
            if(system_switch['ccsd(t)'] == 1): S = S[:-1]
            else:
                print "\n\n            This should not happen.\n"
                printterm(term) 
                set_trace()
        S += '}]'
    return S
        


                  

#********************************
#     Spin-adapted Equations    *
#********************************


                             
""" the weird thing is that this function should not be called in a loop; instead it should be called one a time!!!!"""



def TCE_eq_input_spin_adapted(terms1, filename, leading_mat_name, leading_mat): 
    """ Output spin-adapted equations in TCE format. """

    terms = TCE_terms_change_mat(terms1, 'implicit')
    S = ""
    inds_h = find_index_s_type(terms, 'h')
    inds_p = find_index_s_type(terms, 'p')
    inds_m = find_index_s_type(terms, 'm')
    inds_e = find_index_s_type(terms, 'e')
    S += "range O = 10;\n"
    S += "range V = 100;\n"
    if(len(inds_m) > 0):
        S += "range M = 10;\n"
    if(len(inds_e) > 0):
        S += "range E = 10;\n"
    S += "\n"
    debug = 1
    if(debug and 0):
        set_trace()
    if inds_h:  # produce the header which decalares the indices
        for num_h, indh in enumerate(inds_h):
            if num_h % 8 == 0:
                S += "index  "
                if(num_h == len(inds_h) - 1):
                    S += indh.type + indh.num + " : O;\n"
                else:
                    S += indh.type + indh.num + ", "
            elif num_h % 8  == 7:
                S += indh.type + indh.num + " : O;\n"
            elif num_h == len(inds_h) - 1:
                S += indh.type + indh.num + " : O;\n"
            else:
                S += indh.type + indh.num + ", "
    if inds_p:
        for num_p in range(len(inds_p)):
            if num_p % 8 == 0:
                S += "index  "
                if num_p == len(inds_p) - 1:
                    S += inds_p[num_p].type + inds_p[num_p].num + " : V;\n"
                else:
                    S += inds_p[num_p].type + inds_p[num_p].num + ", "
            elif(num_p % 8 == 7):
                S += inds_p[num_p].type + inds_p[num_p].num + " : V;\n"
            elif(num_p == len(inds_p) - 1):
                S += inds_p[num_p].type + inds_p[num_p].num + " : V;\n"
            else:
                S += inds_p[num_p].type + inds_p[num_p].num + ", "
    if inds_m:
        for num_m in range(len(inds_m)):
            if(num_m % 8 == 0):
                S += "index  "
                if(num_m == len(inds_m) - 1):
                    S += inds_m[num_m].type + inds_m[num_m].num + " : M;\n"
                else:
                    S += inds_m[num_m].type + inds_m[num_m].num + ", "
            elif(num_m % 8 == 7):
                S += inds_m[num_m].type + inds_m[num_m].num + " : M;\n"
            elif(num_m == len(inds_m) - 1):
                S += inds_m[num_m].type + inds_m[num_m].num + " : M;\n"
            else:
                S += inds_m[num_m].type + inds_m[num_m].num + ", "
    if inds_e:
        for num_e in range(len(inds_e)):
            if num_e % 8 == 0:
                S += "index  "
                if(num_e == len(inds_e) - 1):
                    S += inds_e[num_e].type + inds_e[num_e].num + " : E;\n"
                else:
                    S += inds_e[num_e].type + inds_e[num_e].num + ", "
            elif num_e % 8 == 7:
                S += inds_e[num_e].type + inds_e[num_e].num + " : E;\n"
            elif num_e == len(inds_e) - 1:
                S += inds_e[num_e].type + inds_e[num_e].num + " : E;\n"
            else:
                S += inds_e[num_e].type + inds_e[num_e].num + ", "
    S += "\n"
    S += "procedure ncc("
    names = TCE_get_var_names(terms) + [leading_mat_name]
    firstline = 'yes'
    for num_name in range(len(names)):
        if num_name % 4 == 0:
            if firstline == 'yes':
                if num_name == len(names) -1:
                    S += 'out ' + names[num_name] + ')=\n'
                    firstline = 'no'
                else:
                    S += 'in ' + names[num_name] + ', '
                    firstline = 'no'
            else:
                if num_name == len(names) -1:
                    S += '              ' + 'out ' + names[num_name] + ')=\n'
                else:
                    S += '              ' + 'in ' + names[num_name] + ', ' 
        else:
            if num_name == len(names) -1:
                S += '\n              ' + 'out ' + names[num_name] + ')=\n'
            elif num_name % 4 == 3:
                S += 'in ' + names[num_name] + ',\n'
            else:
                S += 'in ' + names[num_name] + ', '
    S += 'begin\n\n'
    S += leading_mat + '==\n'
    for XY in range(len(terms)):
        one_term = TCE_output_one_term(terms[XY])
        if XY == 0 and one_term[0] == '+':
            S += one_term[1:]          # eliminate the '+' sigh for the first term due to some Ohio rule
        else:
            S += one_term
        if XY == len(terms) -1 :
            S += ';\n\n'
        else:
            S += '\n'
    if(0 and system_switch['ccsd(t)']):  # we switch it off                 
        S += 'E_vvvooo[p1,p2,p3,h1,h2,h3]==\n'
        S += 'E_vvvooo[p1,p2,p3,h1,h2,h3]/D[p1,p2,p3,h1,h2,h3]\n\n'
    S += "end"
    if(0 and system_switch['ccsd(t)']):  # we switch it off
        S += '\n\nbegin\n\n'
        S += 'E_tri==\n'
        S += '+1.0*sum[E_vvvooo[p1,p2,p3,h1,h2,h3],{p1,p2,p3,h1,h2,h3}];\n\n'
        S += 'end' 
    file = open(filename, 'w')
    file.write(S)
    file.close()
    return



            

#********************************
#     Spin Orbital Equations    *
#********************************




                             
""" the weird thing is that this function should not be called in a loop; instead it should be called one a time!!!!"""



def TCE_eq_input_spin_orbital(terms1, filename, leading_mat_name, leading_mat):    
    """ Output spin-orbital equations in TCE format. """

    terms = TCE_terms_change_mat(terms1) 
    S = ""
    inds_h = find_index_s_type(terms, 'h')
    inds_p = find_index_s_type(terms, 'p')
    inds_m = find_index_s_type(terms, 'm')
    inds_e = find_index_s_type(terms, 'e')
    S += "range O = 10;\n"
    S += "range V = 100;\n"
    if inds_m: S += "range M = 10;\n"
    if inds_e: S += "range E = 10;\n"
    S += "\n"
    debug = 0
    if debug and 0: set_trace()
    if inds_h:
        for num_h, indh in enumerate(inds_h):
            if num_h % 8  == 0:
                S += "index  "
                if num_h == len(inds_h) - 1: S += indh.type + indh.num + " : O;\n"
                else: S += indh.type + indh.num + ", "
            elif num_h % 8 == 7:
                S += indh.type + indh.num + " : O;\n"
            elif num_h == len(inds_h) - 1:
                S += indh.type + indh.num + " : O;\n"
            else:
                S += indh.type + indh.num + ", "
    if inds_p:
        for num_p, indp in enumerate(inds_p):
            if num_p % 8 == 0:
                S += "index  "
                if num_p == len(inds_p) - 1:
                    S += indp.type + indp.num + " : V;\n"
                else:
                    S += indp.type + indp.num + ", "
            elif num_p % 8 == 7:
                S += indp.type + indp.num + " : V;\n"
            elif num_p == len(inds_p) - 1:
                S += indp.type + indp.num + " : V;\n"
            else:
                S += indp.type + indp.num + ", "
    if inds_m:
        for num_m, indm in enumerate(inds_m):
            if num_m % 8  == 0:
                S += "index  "
                if num_m == len(inds_m) - 1:
                    S += indm.type + indm.num + " : M;\n"
                else:
                    S += indm.type + indm.num + ", "
            elif(num_m % 8 == 7):
                S += indm.type + indm.num + " : M;\n"
            elif(num_m == len(inds_m) - 1):
                S += indm.type + indm.num + " : M;\n"
            else:
                S += indm.type + indm.num + ", "
    if inds_e:
        for num_e, inde in enumerate(inds_e):
            if num_e % 8  == 0:
                S += "index  "
                if num_e == len(inds_e) - 1:
                    S += inde.type + inde.num + " : E;\n"
                else:
                    S += inde.type + inde.num + ", "
            elif num_e % 8  == 7:
                S += inde.type + inde.num + " : E;\n"
            elif num_e == len(inds_e) - 1:
                S += inde.type + inde.num + " : E;\n"
            else:
                S += inde.type + inde.num + ", "
    S += "\n"
    names = TCE_get_var_names(terms) + [leading_mat_name]
    #S += anti_sym(names1)
    S += anti_sym(names)
    S += "procedure ncc("
    firstline = 'yes'
    for num_name, name in enumerate(names):
        if num_name % 4 == 0:
            if firstline == 'yes':
                if num_name == len(names) -1:
                    S += 'out ' + name + ')=\n'
                    firstline = 'no'
                else:
                    S += 'in ' + name + ', '
                    firstline = 'no'
            else:
                if num_name == len(names) -1: S += '              ' + 'out ' + name + ')=\n'
                else: S += '              ' + 'in ' + name + ', ' 
        else:
            if num_name == len(names) -1: S += '\n              ' + 'out ' + name + ')=\n'
            elif num_name % 4  == 3: S += 'in ' + name + ',\n'
            else: S += 'in ' + name + ', '
    S += 'begin\n\n'
    S += leading_mat + '==\n'
    for XY in range(len(terms)):
        S += TCE_output_one_term(terms[XY])
        if XY == len(terms) -1: S += ';\n\n'
        else: S += '\n'
    S += "end"
    file = open(filename, 'w')
    file.write(S)
    file.close()
    return




def anti_sym(names1):                    # claim the antisymmetry
    S = ""
    for xx in names1: S += anti_sym_help(xx[1: xx.index('_')], xx[xx.index('_') + 1: xx.index('[')])
    return S



def anti_sym_help(spin, PH):
    if spin == '':   # implicit spin; in principle, we can simply modify the code below to do this, but let's do it when needed.
        a = b = 0
        if PH[0] == PH[1]: a = 1   # both particle or hole
        if PH[2] == PH[3]: b = 1
        if a != 1 and b != 1: return ""
        else:
            S = "symmetry sym_" + spin + '_' + PH + "[i1:" + PH[0] + ",i2:" + PH[1] + ",i3:" + PH[2] + ",i4:" + PH[3] + "] =\n begin\n"
            if a == 1: S += "  antisym i1; i2\n"
            if b == 1: S += "  antisym i3; i4\n"
            return S + " end\n\n"
    elif spin in ['ab', 'ba', 'a', 'b']: return ""
    elif spin in ['aaa', 'bbb', 'aab', 'abb', 'aa', 'bb']:
        if len(spin) == 3:
            uplist = [PH[0], PH[1], PH[2]]
            lowlist = [PH[3], PH[4], PH[5]]
        if len(spin) == 2:
            uplist = [PH[0], PH[1]]
            lowlist = [PH[2], PH[3]]
        up_ao = []
        up_av = []
        up_bo = []
        up_bv = []
        up_ind_dic = {'ao': up_ao, 'av': up_av, 'bo': up_bo, 'bv': up_bv}
        for xx, yy in enumerate(uplist):
            up_ind_dic[spin[xx] + yy] += [str(xx + 1)]
        low_ao = []
        low_av = []
        low_bo = []
        low_bv = []
        low_ind_dic = {'ao': low_ao, 'av': low_av, 'bo': low_bo, 'bv': low_bv}
        for xx, yy in enumerate(lowlist):
            low_ind_dic[spin[xx] + yy] += [str(xx + 1 + len(spin))]
        if len(up_ao) > 1 or len(up_av) > 1 or len(up_bo) > 1 or len(up_bv) > 1 or len(low_ao) > 1 or len(low_av) > 1 or len(low_bo) > 1 or len(low_bv) > 1:
            if len(spin) == 3:
                S = "symmetry sym_" + spin + '_' + PH + "[i1:" + PH[0] + spin[0] + ",i2:" + PH[1] + spin[1] + ",i3:" + PH[2] + spin[2] + \
                                                        ",i4:" + PH[3] + spin[0] + ",i5:" + PH[4] + spin[1] +  ",i6:" + PH[5] + spin[2] +\
                                                        "] =\n begin\n"
            elif len(spin) == 2:
                S = "symmetry sym_" + spin + '_' + PH + "[i1:" + PH[0] + spin[0] + ",i2:" + PH[1] + spin[1] + \
                                                        ",i3:" + PH[2] + spin[0] + ",i4:" + PH[3] + spin[1] +\
                                                        "] =\n begin\n"
            for uu in [up_ao, up_av, up_bo, up_bv]:
                #if spin == 'abb': set_trace()
                if len(uu) > 1:
                     S += "  antisym "
                     for xx in uu: S += 'i' + xx + '; '
                     S = S[:-2] + "\n"
            for uu in [low_ao, low_av, low_bo, low_bv]:
                if len(uu) > 1:
                     S += "  antisym "
                     for xx in uu: S += 'i' + xx + '; '
                     S = S[:-2] + "\n"
            return S + " end\n\n"         
        else: return ""
    else:
         print spin, PH
         set_trace()
         return




            
#********************************
#     Spin Explicit Equations   *
#********************************


""" the weird thing is that this function should not be called in a loop; instead it should be called one a time!!!!"""


def TCE_eq_input_spin_exp(terms1, filename, leading_mat_name, leading_mat):              
    terms = TCE_terms_change_mat(terms1, 'explicit')     
    S = ""
    inds_ha = find_index_s_type(terms, 'H')
    inds_hb = find_index_s_type(terms, 'h') 
    inds_pa = find_index_s_type(terms, 'P')
    inds_pb = find_index_s_type(terms, 'p') 
    inds_ma = find_index_s_type(terms, 'M')
    inds_mb = find_index_s_type(terms, 'm')
    inds_ea = find_index_s_type(terms, 'E')
    inds_eb = find_index_s_type(terms, 'e')
    S += "range Oa = 10;\n"
    S += "range Ob = 10;\n"
    S += "range Va = 100;\n"
    S += "range Vb = 100;\n"
    if len(inds_ma) > 0: S += "range Ma = 10;\n"
    if len(inds_mb) > 0: S += "range Mb = 10;\n"
    if len(inds_ea) > 0: S += "range Ea = 10;\n"
    if len(inds_eb) > 0: S += "range Eb = 10;\n"
    S += "\n"
    if inds_ha:
        for num_h, indha in enumerate(inds_ha):
            if num_h % 8 == 0:
                S += "index  "
                if num_h == len(inds_ha) - 1:
                    S += transform_exp_spin_ind_to_tce_format(indha) + " : Oa;\n"                      # transform indices from  'capital-small' to 'ab' notation  
                else:
                    S += transform_exp_spin_ind_to_tce_format(inds_ha[num_h]) + ", "
            elif num_h % 8 == 7:
                S += transform_exp_spin_ind_to_tce_format(indha) + " : Oa;\n"
            elif num_h == len(inds_ha) - 1:
                S += transform_exp_spin_ind_to_tce_format(indha) + " : Oa;\n"
            else:
                S += transform_exp_spin_ind_to_tce_format(indha) + ", "
    if inds_hb:
        for num_h, indhb in enumerate(inds_hb):
            if num_h % 8 == 0:
                S += "index  "
                if num_h == len(inds_hb) - 1:
                    S += transform_exp_spin_ind_to_tce_format(indhb) + " : Ob;\n"                      # transform indices from  'capital-small' to 'ab' notation
                else:
                    S += transform_exp_spin_ind_to_tce_format(indhb) + ", "
            elif(num_h % 8 == 7):
                S += transform_exp_spin_ind_to_tce_format(indhb) + " : Ob;\n"
            elif num_h == len(inds_hb) - 1:
                S += transform_exp_spin_ind_to_tce_format(indhb) + " : Ob;\n"
            else:
                S += transform_exp_spin_ind_to_tce_format(indhb) + ", "
    if inds_pa:
        for num_p in range(len(inds_pa)):
            if num_p % 8  == 0:
                S += "index  "
                if num_p == len(inds_pa) - 1:
                    S += transform_exp_spin_ind_to_tce_format(inds_pa[num_p]) + " : Va;\n"
                else:
                    S += transform_exp_spin_ind_to_tce_format(inds_pa[num_p]) + ", "
            elif num_p % 8  == 7:
                S += transform_exp_spin_ind_to_tce_format(inds_pa[num_p]) + " : Va;\n"
            elif num_p == len(inds_pa) - 1:
                S += transform_exp_spin_ind_to_tce_format(inds_pa[num_p]) + " : Va;\n"
            else:
                S += transform_exp_spin_ind_to_tce_format(inds_pa[num_p]) + ", "
    if inds_pb:
        for num_p in range(len(inds_pb)):
            if num_p % 8 == 0:
                S += "index  "
                if num_p == len(inds_pb) - 1:
                    S += transform_exp_spin_ind_to_tce_format(inds_pb[num_p]) + " : Vb;\n"
                else:
                    S += transform_exp_spin_ind_to_tce_format(inds_pb[num_p]) + ", "
            elif num_p % 8 == 7:
                S += transform_exp_spin_ind_to_tce_format(inds_pb[num_p]) + " : Vb;\n"
            elif num_p == len(inds_pb) - 1:
                S += transform_exp_spin_ind_to_tce_format(inds_pb[num_p]) + " : Vb;\n"
            else:
                S += transform_exp_spin_ind_to_tce_format(inds_pb[num_p]) + ", "
    if inds_ma:
        for num_m in range(len(inds_ma)):
            if num_m % 8 == 0:
                S += "index  "
                if(num_m == len(inds_ma) - 1):
                    S += transform_exp_spin_ind_to_tce_format(inds_ma[num_m]) + " : Ma;\n"
                else:
                    S += transform_exp_spin_ind_to_tce_format(inds_ma[num_m]) + ", "
            elif num_m % 8 == 7:
                S += transform_exp_spin_ind_to_tce_format(inds_ma[num_m]) + " : Ma;\n"
            elif num_m == len(inds_ma) - 1:
                S += transform_exp_spin_ind_to_tce_format(inds_ma[num_m]) + " : Ma;\n"
            else:
                S += transform_exp_spin_ind_to_tce_format(inds_ma[num_m]) + ", "
    if inds_mb:
        for num_m in range(len(inds_mb)):
            if num_m % 8 == 0:
                S += "index  "
                if num_m == len(inds_mb) - 1:
                    S += transform_exp_spin_ind_to_tce_format(inds_mb[num_m]) + " : Mb;\n"
                else:
                    S += transform_exp_spin_ind_to_tce_format(inds_mb[num_m]) + ", "
            elif num_m % 8 == 7:
                S += transform_exp_spin_ind_to_tce_format(inds_mb[num_m]) + " : Mb;\n"
            elif num_m == len(inds_mb) - 1:
                S += transform_exp_spin_ind_to_tce_format(inds_mb[num_m]) + " : Mb;\n"
            else:
                S += transform_exp_spin_ind_to_tce_format(inds_mb[num_m]) + ", "
    if inds_ea:
        for num_e in range(len(inds_ea)):
            if num_e % 8 == 0:
                S += "index  "
                if(num_e == len(inds_ea) - 1):
                    S += transform_exp_spin_ind_to_tce_format(inds_ea[num_e]) + " : Ea;\n"
                else:
                    S += transform_exp_spin_ind_to_tce_format(inds_ea[num_e]) + ", "
            elif num_e % 8 == 7:
                S += transform_exp_spin_ind_to_tce_format(inds_ea[num_e]) + " : Ea;\n"
            elif num_e == len(inds_ea) - 1:
                S += transform_exp_spin_ind_to_tce_format(inds_ea[num_e]) + " : Ea;\n"
            else:
                S += transform_exp_spin_ind_to_tce_format(inds_ea[num_e]) + ", "
    if inds_eb:
        for num_e in range(len(inds_eb)):
            if num_e % 8  == 0:
                S += "index  "
                if(num_e == len(inds_eb) - 1):
                    S += transform_exp_spin_ind_to_tce_format(inds_eb[num_e]) + " : Eb;\n"
                else:
                    S += transform_exp_spin_ind_to_tce_format(inds_eb[num_e]) + ", "
            elif num_e % 8 == 7:
                S += transform_exp_spin_ind_to_tce_format(inds_eb[num_e]) + " : Eb;\n"
            elif num_e == len(inds_eb) - 1:
                S += transform_exp_spin_ind_to_tce_format(inds_eb[num_e]) + " : Eb;\n"
            else:
                S += transform_exp_spin_ind_to_tce_format(inds_eb[num_e]) + ", "
    S += "\n"
    names1 = TCE_get_var_names(terms)              # give the mat names in spin-explicit form                       
    names = names1 + [leading_mat_name]
    S += anti_sym(names1)                          # claim the antisymmetry
    S += "procedure ncc("
    firstline = 'yes'
    for num_name in range(len(names)):             # tells the input and output quantity
        if num_name % 4 == 0:
            if firstline == 'yes':
                if num_name == len(names) -1:
                    S += 'out ' + names[num_name] + ')=\n'
                    firstline = 'no'
                else:
                    S += 'in ' + names[num_name] + ', '
                    firstline = 'no'
            else:
                if num_name == len(names) -1:
                    S += '              ' + 'out ' + names[num_name] + ')=\n'
                else:
                    S += '              ' + 'in ' + names[num_name] + ', ' 
        else:
            if num_name == len(names) -1:
                S += '\n              ' + 'out ' + names[num_name] + ')=\n'
            elif num_name % 4 == 3:
                S += 'in ' + names[num_name] + ',\n'
            else:
                S += 'in ' + names[num_name] + ', '
    S += 'begin\n\n'
    S += leading_mat + '==\n'
    for XY in range(len(terms)):
        S += TCE_output_one_term_spin_exp(terms[XY])        # output terms  
        if XY == len(terms) -1: S += ';\n\n'
        else: S += '\n'
    S += "end"
    file = open(filename, 'w')
    file.write(S)
    file.close()
    return





                                                                                                                  
def transform_exp_spin_ind_to_tce_format(ind):
    dic_1 = {'P': 'a',
             'H': 'a',
             'M': 'a',
             'E': 'a',
             'p': 'b',
             'h': 'b',
             'm': 'b',
             'e': 'b'}
    return PHME_2_phme(ind.type) + ind.num + dic_1[ind.type] 


    

def TCE_output_one_term_spin_exp(term):                                                 # output one term A*B*C in spin-explicit form
    S = ""
    inds = []
    for xx in term.coefficient.matElement:
        inds += xx.matUpperIndicees + xx.matLowerIndicees
    contracted_inds = find_contracted_inds(inds)
    contracted_inds = reorder_inds_hmpe(contracted_inds)
    if term.coefficient.const > 0: S += '+' + str(term.coefficient.const)
    elif term.coefficient.const < 0: S += str(term.coefficient.const)
    else:
        print "\n\n                                       Zero coefficient? \n"
        set_trace()
    if len(term.coefficient.matElement) == 0:
        print "No matelement?".center(60)
        set_trace()
    elif len(contracted_inds) == 0:
        for mat in term.coefficient.matElement:
            S += "*" + mat.name + '['
            clock = 0
            for indd in mat.matUpperIndicees + mat.matLowerIndicees:
                if clock == len(mat.matUpperIndicees + mat.matLowerIndicees) -1:
                    S += transform_exp_spin_ind_to_tce_format(indd)    
                else:
                    S += transform_exp_spin_ind_to_tce_format(indd)+ ','
                clock += 1
            S += ']'
    else:
        S += '*sum['
        clock1 = 0
        for mat in term.coefficient.matElement:
            S += mat.name + '['
            clock2 = 0
            for indd in mat.matUpperIndicees + mat.matLowerIndicees:
                if clock2 == len(mat.matUpperIndicees + mat.matLowerIndicees) -1:
                    S += transform_exp_spin_ind_to_tce_format(indd)
                else:
                    S += transform_exp_spin_ind_to_tce_format(indd)+ ','
                clock2 += 1
            S += ']'
            if clock1 == len(term.coefficient.matElement)-1: S += ','
            else: S += '*'
            clock1 += 1
        S += '{'
        clock3 = 0
        for contracted_ind in contracted_inds:
            S += transform_exp_spin_ind_to_tce_format(contracted_ind)
            if clock3 != len(contracted_inds)-1:
                S += ',' 
            clock3 += 1
        S += '}]'
    return S

