#! /usr/bin/python
""" Derives equations for SF-[2]_R12 correction in J. Chem. Phys. 135, 214105 (2011).
    
   One keypoint to remember: during the derivation, we only used two orb sets: hole and virtual orbitals.
   Hole is to simulate OBS; virtual is to simulate CABS; the constraint of correlating only certain orbs (subset of occupied orbs)
   is automatically taken care of by the geminal amplitudes, which are nonzero only in predefined GG/gg space.
   I get confused sometimes. The derivation is correct anyway """

from MBEq import *
from pdb import set_trace
import copy, pickle, sys, os


screenwidth = 90
os.system('date')
print "\n\n"
#set_trace()


system_switch['normal order'] = 'V'
system_switch['use spin orbital'] = 0 
system_switch['explicit spin'] = 0
system_switch['GWT CU truncate rank'] = 1
system_switch['Allow Internal Contraction'] = 0
system_switch['canonical scheme'] = 'p'


# V 

if 1:
    [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12] = hfactory(12)
    [pp1, pp2, pp3] = pfactory(3)
    [g1, g2, g3] = gfactory(3)
    [r, s, p, q, v, w, t, u, x, y, z] = [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11]
    alpha_p = pp1
    beta_p = pp2
    alpha2_p = pp2
    alpha3_p = pp3
    p2 = h9
    p3 = h10
    q2 = h11
    q3 = h12
    [kappa, lambdaa, ksi] = [g1, g2, g3]
    R_1 = direct_term_construction([0.5, [ ['t2', [p, q, r, s]], ['r', [r, s, alpha_p, beta_p]] ],[ [alpha_p, beta_p, p, q] ] ])
    R_2 = direct_term_construction([1.0, [ ['t2', [p, q, r, s]], ['r', [r, s, alpha_p, x]] ],[ [alpha_p, x, p, q] ] ])
    R_3 = direct_term_construction([-1.0, [ ['t2', [p, q, r, s]], ['r', [r, s, alpha_p, x]], ['Gamma', [z, x, p, q]], ['x1', [y, z]] ],[ [alpha_p, y] ] ])
    R = [R_1, R_2, R_3]
    #printterms([R_1, R_2, R_3])

#*********************   V TERM **********************
    # V term
    if 1:  # 2 <0| [H, \Omega^1] | 0 >              
        H2 = H()[1:]                    # it has been verified that the core H does not contribute
        H2 = expand_general_indices(H2)
        H2 = merge_all(H2)
        V = commutator(H2, R)      # 'eta' has been replaced by delta here to emulate genuine vacuum NO
        V = multiply_co(2.0, V)   # since it is 2 * ....
        V = mc_density(V)
        V = replace_sf_3_dm_with_cumu_s(V)
        V = terms_changemat_s_name(V, ['d2', 'd1', 'Gamma1'], ['Gamma', 'Gamma1', 'Lambda1'])
        V = merge_all(V)
        V = replace_2b_dm_with_cumu_s(V)
        V = mr_f12_valeev_discard_quadra_cumu(V)
        V = merge_all(V)               
        V = combine_lambda1_and_inverse_to_delta_s(V)
        V = km_gwt_decode_delta_s(V)
        V = merge_all(V)
        V = mr_f12_valeev_discard_s(V)

        V = replace_2b_cu_with_specified_inds_terms(V, [])   # p q
        V = merge_all(V)              
        #set_trace()
        V = valeev_f12_resolve_2cabs_terms(V) 
        #set_trace()
        V = merge_all(V)                # 3 terms
        common_factors = V[0].coefficient.matElement[:1]
        V = separate_terms_wrt_factorization(V, common_factors)[0][1]
        V = group_terms_based_on_matname_list(V, ['Gamma'])
        V = group_factorize_terms_based_on_mats_list(V, [V[0].coefficient.matElement[:1]])
        termslist_texform('sf', V, 'Vfac.tex', 'V^{rs}_{pq}', "", 'sf:valeev:mrf12:new', 3, 1, 0)
  


#*********************   B TERM  **********************

# since the maximum RDM rank is 3
if 1:
    [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12] = hfactory(12)
    [pp1, pp2, pp3] = pfactory(3)
    [g1, g2, g3] = gfactory(3)
    [r, s, p, q, v, w, t, u, x, y, z] = [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11]
    if 0:
        for x in [h1, h2, h3, h4, h5, h6, h7, h8]: x.fix = 1
    alpha_p = pp1
    beta_p = pp2
    alpha2_p = pp2
    alpha3_p = pp3
    p2 = h9
    p3 = h10 
    q2 = h11
    q3 = h12 
    [kappa, lambdaa, ksi] = [g1, g2, g3] 
    R_1 = direct_term_construction([0.5, [ ['t2', [p, q, r, s]], ['r', [r, s, alpha_p, beta_p]] ],[ [alpha_p, beta_p, p, q] ] ])
    R_2 = direct_term_construction([1.0, [ ['t2', [p, q, r, s]], ['r', [r, s, alpha_p, x]] ],[ [alpha_p, x, p, q] ] ])
    R_3 = direct_term_construction([-1.0, [ ['t2', [p, q, r, s]], ['r', [r, s, alpha_p, x]], ['Gamma', [z, x, p, q]], ['x1', [y, z]] ],[ [alpha_p, y] ] ])
    R = [R_1, R_2, R_3]
    # since R/L only contain 'Lambda' without 'Lambda1', so they won't be discarded in the tricky step of solvetwoop()
    L_1 = direct_term_construction([0.5, [ ['t2', [v, w, t, u]], ['r', [alpha_p, beta_p, v, w]] ],[ [t, u, alpha_p, beta_p] ] ])
    L_2 = direct_term_construction([1.0, [ ['t2', [v, w, t, u]], ['r', [alpha_p, x, v, w]] ],[ [t, u, alpha_p, x] ] ])
    L_3 = direct_term_construction([-1.0, [ ['t2', [v, w, t, u]], ['r', [alpha_p, x, v, w]], ['Gamma', [t, u, z, x]], ['x1', [z, y]] ],[ [y, alpha_p] ] ])
    L = [L_1, L_2, L_3]
    Lpp = L[:1]
    Lph = L[1:]
    Rpp = R[:1]
    Rph = R[1:]
    #printterms([R_1, R_2, R_3, L_1, L_2, L_3])

    F = H()[:1]                    # Fock
    F = expand_general_indices(F)


    # (1,1) + 2 (1, 23) + (23, 23)
    if 1:                              
        X = commutator(L[:], commutator(F, R[:]))        
        X = mc_density(X)             
        X = merge_all(X)              
        X = replace_sf_3_dm_with_cumu_s(X)  
        X = replace_2b_cumu_with_dm_s(X)  
        X = terms_changemat_s_name(X, ['d3', 'd2', 'd1', 'Gamma1'], ['Gamma', 'Gamma', 'Gamma1', 'Lambda1'])
        X = merge_all(X)   
        X = replace_2b_dm_with_cumu_s(X) 
        X = merge_all(X)   
        X = combine_lambda1_and_inverse_to_delta_s(X)
        X = km_gwt_decode_delta_s(X)
        X = mr_f12_valeev_discard_s(X) 
        X = merge_all(X)               # 37
        Xrdm = replace_2b_cu_with_specified_inds_terms(X, [])   # since we know that there will be just 2-rdm in this case, we simply replace all cumu by rdms. 
        Xrdm = merge_all(Xrdm)         # 21
        Xrdm = valeev_f12_resolve_2cabs_terms(Xrdm) # 98
        Xrdm = merge_all(Xrdm)         # 22           
        print len(Xrdm)
        texform('sf', Xrdm, 'Brdm.tex', '', "", 'sf:valeev:mrf12:new', 3)
        xrdmfile = open('Xrdm', 'w')
        pickle.dump(Xrdm, xrdmfile)
        xrdmfile.close()





    if 1:
        xrdmfile = open('Xrdm', 'r')
        X = pickle.load(xrdmfile)
        sep_tt =  2
        if 1:
            ts = X[0].coefficient.matElement[:2]
            X = separate_terms_wrt_factorization(X, ts)[0][1]
            sep_tt = 0
        X = group_terms_based_on_matname_list(X, ['Gamma'])
        #texform('sf', X, 'Bgamma.tex', '', "", 'sf:valeev:mrf12:new', 3)
        #set_trace()
        #Xfac = group_factorize_terms_based_on_mats_list(X, [X[0].coefficient.matElement[:(1+sep_tt)], X[6].coefficient.matElement[:(2+ sep_tt)], X[9].coefficient.matElement[:(1+sep_tt)], X[0].coefficient.matElement[:sep_tt]])
        Xfac = group_factorize_terms_based_on_mats_list(X, [X[0].coefficient.matElement[:(1+sep_tt)], X[4].coefficient.matElement[:(2+ sep_tt)], X[7].coefficient.matElement[:(3+ sep_tt)]])
        termslist_texform('sf', Xfac, 'Bfac.tex')


    if 1:
        xrdmfile = open('Xrdm', 'r')
        X = pickle.load(xrdmfile)
        X = group_terms_based_on_matname_list(X, ['Gamma'])
        printterms(X)
        sep_tt =  2
        mats1 = X[0].coefficient.matElement[:(1+ sep_tt)] 
        mats2 = X[4].coefficient.matElement[:(2+ sep_tt)]
        mats3 = [X[9].coefficient.matElement[0], X[9].coefficient.matElement[5]] 
        mats4 = [X[7].coefficient.matElement[0], X[7].coefficient.matElement[5]] 
        facmatlist = [mats1, mats2, mats3, mats4]
        Xfac = group_factorize_terms_based_on_mats_list(X, facmatlist)
        termslist_texform('sf', Xfac, 'Bfac2.tex', 'H^{(2)}')







