'''
        $$$ Cálculo das condições psicométricas do ar úmido $$$$
Variáveis importantes:
   patm  - pressão barométrica local (kPa)
   tbs   - Temperatura do bulbo seco (°C)
   ur    - Umidade relativa (%)
   pvs   - Pressão de vapor de saturação (kPa)
   pv    - Pressão parcial de vapor (kPa)
   rm    - Razão de mistura (%)
   tpo   - Temperatura do ponto de orvalho (°C)
   tbm   - Temperatura do bulbo molhado (°C)
   e     - Enthalpy (kJ/kg)
   ve    - Volume específico (m3/kg)
   q     - Vazão total de ar  (m3/hora)
   
   Não há validação dos dados de entrada
'''

def pressao_vapor_saturado(t):
    # cálculo da pressão do vapor de saturação
    # log in Python é ln na matemática
    t = t + 273.16
    if t > 273.16:
        aux = -7511.52 / t + 89.63121 + 0.023998970 * t
        aux = aux - 1.1654551E-5 * (t ** 2) - 1.2810336E-8 * (t ** 3)
        aux = aux + 2.0998405E-11 * (t ** 4) - 12.150799 * np.log(t)
        p_vs = np.exp(aux)
        return p_vs
    else:
        aux = 24.2779 - 6238.64 / t - 0.344438 * np.log(t)
        p_vs = np.exp(aux)
        return p_vs

def razao_mistura1(p):
    # Primeiro método de cálculo para razão de mistura
    global patm
    rm_1 = 0.62198 * p / (patm - p)
    return rm_1

def razao_mistura2(ts, th, w):
    # Segundo  método de cálculo para razão de mistura
    aux1 = (2501. - 2.411 * th) * w - 1.006 * (ts - th)
    aux2 = 2501. + 1.775 * ts - 4.186 * th
    rm_2 = aux1 / aux2
    return rm_2

def umidade_relativa(p, ps):
    # Cálculo da Umidade Relativa
    u_rel = p / ps
    return u_rel

def entalpia(t, w):
    # Cálculo de Entalpia
    h = 1.006*t + w*(2501.+1.775*t)
    return h

def pressao_vapor(w):
    # Cálculo da pressão parcial de vapor
    global patm
    p_vap = patm * w / (0.62198 + w)
    return p_vap

def temperatura_ponto_orvalho(p):
    # Cálculo da temperatura do ponto de orvalho
    a = np.log10(p * 10)  # log10 requires Numpy
    t_po = (186.4905 - 237.3 * a) / (a - 8.2859)
    return t_po

def temperatura_b_molhado(ts,et):
    # Cálculo da temperatura do bulbo molhado
    global patm
    delta = 0.1
    th = ts-delta
    while True:
        rmbs=(et-1.006*th)/(2501.+1.775*th)
        ps=pressao_vapor_saturado(th)
        urel=(patm*rmbs)/(ps*(0.62198+rmbs))
        if urel>1:
            th=th+delta
            delta=delta/2
        if urel<0.999:
            th=th-delta
        if 0.999<=urel<1:
            break
    t_bm=th
    return t_bm

def temperatura_b_seco(h,w):
    # Cálculo da temperatura do bulbo seco
    t_bs=(h-2501.*w)/(1.006+1.775*w)
    return t_bs

def volume_especifico(t,w):
    # Cálculo de volume específico
    global patm
    r=0.28705
    t=t+273.16
    v_esp=r*t/patm*(1+1.6078*w)
    return v_esp

def print_hi(name):
    print(f"  {name}")

def p_atm():
    # Cálculo da pressão barométrica, kPa
    alt = float(input('Altitude (m)......................... '))
    a = 2.2556e-5
    b = 5.2559
    pre_atm=101.324*(1-a*alt) ** b
    return pre_atm

def pe_tbs_ur():
    # Ponto de Estado  f (tbs, ur)
    global tbs, tbm, patm,ur, pvs, pv, rm, tpo, ve, e
    tbs = float(input('Temperatura de bulbo seco (ºC)....... '))  
    ur =  float(input('Umidade relativa (%) ................ '))  
    ur = ur / 100.

    if ur == 1.:
        ur = 0.99999
    pvs = pressao_vapor_saturado(tbs)
    pv = ur * pvs
    rm = razao_mistura1(pv)
    e = entalpia(tbs, rm)
    ve = volume_especifico(tbs, rm)
    if ur == 0.99:
        tpo = tbs
    else:
        tpo = temperatura_ponto_orvalho(pv)
    tbm = temperatura_b_molhado(tbs, e)

def pe_tbs_tbm():
    # Ponto de Estado -   f (tbs, tbm)
    global tbs, tbm, patm,ur, pvs, pv, rm, tpo, ve, e
    tbs = float(input('Temperatura de bulbo seco (ºC)....... '))  
    tbm = float(input('Temperatura de bulbo molhado (ºC).... '))  
    pvs=pressao_vapor_saturado(tbs)
    if tbs == tbm:
        ur=1.
        tpo=tbs
        pv=pvs
        rm=razao_mistura1(pvs)
    else:
        pvsu=pressao_vapor_saturado(tbm)
        rmsu=razao_mistura1(pvsu)
        rm= razao_mistura2(tbs,tbm,rmsu)
        pv=pressao_vapor(rm)
        ur=umidade_relativa(pv,pvs)
        tpo=temperatura_ponto_orvalho(pv)
    e=entalpia(tbs,rm)
    ve =volume_especifico(tbs,rm)

def pe_tbs_tpo():
    # Ponto de Estado -   f (tbs, tpo)
    global tbs, tbm, patm, tbs, ur, pvs, pv, rm, tpo, ve, e
    tbs= float(input('Temperatura de bulbo seco (ºC)............. '))  
    tpo= float(input('Temperatura de ponto de orvalho (ºC)....... '))  
    pvs=pressao_vapor_saturado(tbs)
    if tbs == tpo:
        pv=pvs
        ur=0.999999
        rm=razao_mistura1(pvs)
        e=entalpia(tbs,rm)
        tbm=tbs
    else:
        pv=pressao_vapor_saturado(tpo)
        rm=razao_mistura1(pv)
        ur=umidade_relativa(pv,pvs)
        e=entalpia(tbs,rm)
        tbm=temperatura_b_molhado(tbs,e)
    ve=volume_especifico(tbs,rm)

def aquece_resfria():
    # Processo 1 - Aquecimento ou resfriamento
    global tbs, tbm, patm, tbs, ur, pvs, pv, rm, tpo, ve, e
    #
    # Dados Pontos de Estado 1 e 2
    print('********* Ponto de Estado 1 **********')
    tbs1=float(input('Temperatura bulbo seco(ºC)... '))
    ur1=float(input('Umidade relativa(%).......... '))
    print('********* Ponto de Estado 2 **********')
    tbs2 =float(input('Temperatura bulbo seco(ºC)... '))
    #
    # Point State 1
    ur=ur1/100.
    pvs=pressao_vapor_saturado(tbs1)
    pv=ur*pvs
    rm=razao_mistura1(pv)
    tpo=temperatura_ponto_orvalho(pv)
    e=entalpia(tbs1,rm)
    tbm=temperatura_b_molhado(tbs1,e)
    ve=volume_especifico(tbs1,rm)
    tbs=tbs1
    # Ponto de Estado 2
    pvs2=pressao_vapor_saturado(tbs2)
    if tbs2 > tpo:
        rm2=rm
        pv2=pv
        ur2=pv2/pvs2
        tpo2=tpo
        e2=entalpia(tbs2,rm2)
        tbm2=temperatura_b_molhado(tbs2,e2)
    else:
        ur2=1.
        pv2=ur2*pvs2
        rm2=razao_mistura1(pv2)
        tpo2=tbs2
        e2=entalpia(tbs2,rm2)
        tbm2=tbs2
    ve2=volume_especifico(tbs2,rm2)
    # Resultados do Ponto de Estado 1 e 2
    #                      Ponto de Estado 1
    print('\n************** Aquecimento / Resfriamento***************')
    print('************** Ponto de Estado 1 *********************')
    print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs)
    print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm)
    print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo)
    print(' Umidade relativa(%%)..................... %7.3f' % ur)
    print(' Razão de mistura (g/kg) ................ %7.3f' % rm)
    print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv)
    print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs)
    print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
    print(' entalpia(kJ/kg) ........................ %7.3f' % e)
    print(' Volume específico(m3/kg) ............... %7.3f' % ve)
    #                      Ponto de Estado 2
    ur2 = ur2 * 100.
    rm2 = rm2 * 1000.
    print('\n************** Ponto de Estado 2 *********************')
    print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs2)
    print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm2)
    print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo2)
    print(' Umidade relativa(%%)..................... %7.3f' % ur2)
    print(' Razão de mistura (g/kg) ................ %7.3f' % rm2)
    print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv2)
    print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs2)
    print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
    print(' entalpia(kJ/kg) ........................ %7.3f' % e2)
    print(' Volume específico(m3/kg) ............... %7.3f' % ve2)

def u_adiabatica_tbs():
    # Processo 2
    # Umidificação adiabática - f(tbs1, ur1, tbs2)
    global tbs, tbm, patm, tbs, ur, pvs, pv, rm, tpo, ve, e
    #
    # Dados Pontos de Estado 1 e 2
    print('********* Ponto de Estado 1 **********')
    tbs1 = float(input('Temperatura bulbo seco(ºC)... '))
    ur1  = float(input('Umidade relativa(%).......... '))
    print('********* Ponto de Estado 2 **********')
    tbs2 = float(input('Temperatura bulbo seco(ºC)... '))
    #
    # Ponto de Estado 1
    ur = ur1/100.
    tbs=tbs1
    pvs = pressao_vapor_saturado(tbs)
    pv  = ur*pvs
    rm  = razao_mistura1(pv)
    e   = entalpia(tbs,rm)
    ve  = volume_especifico(tbs,rm)
    if ur == 1:
        tpo = tbs
        tbm = tbs
    else:
        tpo = temperatura_ponto_orvalho(pv)
        tbm = temperatura_b_molhado(tbs,e)
    # Ponto de Estado 2
    tbm2 = tbm
    e2   = e
    rm2  = rm
    while True:
        if tbs > tbs2:  # indica que ur1 < ur2
            rm2 = rm2+0.0001
        else:
            rm2 = rm2-0.0001
        pv2  = pressao_vapor(rm2)
        tbs0 = temperatura_b_seco(e2,rm2)
        pvs2 = pressao_vapor_saturado(tbs2)
        ur2  = pv2/pvs2
        if tbs > tbs2:   # tbs > tbs2
            break
        elif tbs0 > tbs2:
            break
    tpo2 = temperatura_ponto_orvalho(pv2)
    ve2  = volume_especifico(tbs2,rm2)
    # Resultados do  Ponto de Estado 1 e 2
    #                      Ponto de Estado 1
    ur = ur * 100.
    rm = rm * 1000.
    print('\n************** Aquecimento / Resfriamento***************')
    print('************** Ponto de Estado 1 *********************')
    print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs)
    print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm)
    print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo)
    print(' Umidade relativa(%%)..................... %7.3f' % ur)
    print(' Razão de mistura (g/kg) ................ %7.3f' % rm)
    print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv)
    print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs)
    print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
    print(' Entalpia(kJ/kg) ........................ %7.3f' % e)
    print(' Volume específico(m3/kg) ............... %7.3f' % ve)
    #                      Ponto de Estado 2
    ur2 = ur2 * 100.
    rm2 = rm2 * 1000.
    print('\n************** Ponto de Estado 2 *********************')
    print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs2)
    print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm2)
    print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo2)
    print(' Umidade relativa(%%)..................... %7.3f' % ur2)
    print(' Razão de mistura (g/kg) ................ %7.3f' % rm2)
    print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv2)
    print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs2)
    print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
    print(' Entalpia(kJ/kg) ........................ %7.3f' % e2)
    print(' Volume específico(m3/kg) ............... %7.3f' % ve2)

def u_adiabatica_ur():
    # Processo 3
    # Umidificação adiabática - f(tbs1, ur1, ur2)
    #
    global tbs, tbm, patm, tbs, ur, pvs, pv, rm, tpo, ve, e
    #
    # Input - Ponto de Estados 1 and 2
    print('********* Ponto de Estado 1 **********')
    tbs = float(input('Temperatura bulbo seco(ºC).... '))
    ur  = float(input('Umidade relativa (%).......... '))
    print('********* Ponto de Estado 2 **********')
    ur2 = float(input('Umidade relativa (%).......... '))
    #
    # Ponto de Estado 1
    ur=ur / 100.
    if ur2 == 1.:
        ur2=0.9999
    pvs= pressao_vapor_saturado(tbs)
    pv=ur * pvs
    rm=razao_mistura1(pv)
    e=entalpia(tbs, rm)
    ve=volume_especifico(tbs,rm)
    tpo=temperatura_ponto_orvalho(pv)
    tbm=temperatura_b_molhado(tbs,e)
    # Ponto de estado 2
    ur2=ur2 / 100.
    e2=e
    tbm2=tbm
    rm2=rm
    while True:
        rm2 = rm2 + 0.0001
        pv2 = pressao_vapor(rm2)
        tbs2 = temperatura_b_seco(E2, rm2)
        pvs2 = pressao_vapor_saturado(tbs2)
        ur0 = pv2 / pvs2
        if ur0 >= ur2:
            break
    if ur2 == 1.:
        tbs2=tbm
        tpo2=tbm
    else:
        tpo2 = temperatura_ponto_orvalho(pv2)
    ve2 = volume_especifico(tbs2,rm2)
    # Resultados do Ponto de Estado 1 e 2
    #                      Ponto de Estado 1
    ur = ur * 100.
    rm = rm * 1000.
    print('\n************** Umidificação Adiabática *****************')
    print('************** Ponto de Estado 1 *********************')
    print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs)
    print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm)
    print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo)
    print(' >Umidade relativa(%%)................... %7.3f' % ur)
    print(' Razão de mistura (g/kg) ................ %7.3f' % rm)
    print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv)
    print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs)
    print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
    print(' Entalpia(kJ/kg) ........................ %7.3f' % e)
    print(' Volume específico(m3/kg) ............... %7.3f' % ve)
    #                      Ponto de Estado 2
    ur2 = ur2 * 100.
    rm2 = rm2 * 1000.
    print('\n************** Ponto de Estado 2 *********************')
    print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs2)
    print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm2)
    print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo2)
    print(' >Umidade relativa(%%).................... %7.3f' % ur2)
    print(' Razão de mistura (g/kg) ................ %7.3f' % rm2)
    print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv2)
    print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs2)
    print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
    print(' Entalpia(kJ/kg) ........................ %7.3f' % e2)
    print(' Volume específico(m3/kg) ............... %7.3f' % ve2)

def u_adiabatica_rm():
    # Processo 4
    # Umidificação adiabática -  f (tbs1, rm1, rm2).
    #
    global tbs, tbm, patm, tbs, ur, pvs, pv, rm, tpo, ve, e
    #
    # Input - Ponto de Estados 1 and 2
    print('********* Ponto de Estado 1 **********')
    tbs=float(input('Temperatura bulbo seco(ºC)... '))
    rm=float(input('Razão de mistura (g/kg)...... '))
    print('********* Ponto de Estado 2 **********')
    rm2=float(input('Razão de mistura (g/kg).......'))
    #
    rm=rm/1000.
    rm2=rm2/1000.
    # Ponto de Estado 1
    pvs= pressao_vapor_saturado(tbs)
    pv=  pressao_vapor(rm)
    ur=(pv/pvs)
    if ur>1.:
        print(' O valor da razão de mistura do ponto 1 é muito alto')
    else:
        e=entalpia(tbs,rm)
        ve=volume_especifico(tbs,rm)
        tpo=temperatura_ponto_orvalho(pv)
        tbm=temperatura_b_molhado(tbs,e)
        # Status Point 2
        e2=e
        tbm2=tbm
        tbs2=temperatura_b_seco(e2,rm2)
        pvs2=pressao_vapor_saturado(tbs2)
        pv2=pressao_vapor(rm2)
        ur2=(pv2/pvs2)
        if ur2>1.:
            print(' O valor da razão de mistura fornecida no ponto 2 é muito alta')
        else:
            tpo2=temperatura_ponto_orvalho(pv2)
            ve2=volume_especifico(tbs2,rm2)
            if ur2 >= 1.0:
                ur2= 0.99999
                tbs2=tbm2
                tpo2=tbm2
                pvs2=pressao_vapor_saturado(tbs2)
                pv2=pvs2
                rm2=razao_mistura1(pv2)
                ve2=volume_especifico(tbs2,rm2)
            # Resultados do Ponto de Estado 1 e 2
            #                      Ponto de Estado 1
            ur=ur*100.
            rm=rm*1000.
            print('\n************** Umidificação Adiabática *****************')
            print('************** Ponto de Estado 1 *********************')
            print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs)
            print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm)
            print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo)
            print(' Umidade relativa(%%)..................... %7.3f' % ur)
            print(' >Razão de mistura (g/kg) ............... %7.3f' % rm)
            print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv)
            print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs)
            print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
            print(' Entalpia(kJ/kg) ........................ %7.3f' % e)
            print(' Volume específico(m3/kg) ............... %7.3f' % ve)
            #                      Ponto de Estado 2
            ur2 = ur2 * 100.
            rm2 = rm2 * 1000.
            print('\n************** Ponto de Estado 2 *********************')
            print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs2)
            print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm2)
            print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo2)
            print(' Umidade relativa(%%)..................... %7.3f' % ur2)
            print(' >Razão de mistura (g/kg) ............... %7.3f' % rm2)
            print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv2)
            print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs2)
            print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
            print(' Entalpia(kJ/kg) ........................ %7.3f' % e2)
            print(' Volume específico(m3/kg) ............... %7.3f' % ve2)

def mistura_fluxos():
    # Processo 5
    # Mistura de dois fluxos de ar
    # Input - Ponto de Estados 1 e 2
    print('********* Ponto de Estado 1 **********')
    tbs1=float(input('Temperatura bulbo seco(ºC)... '))
    ur1= float(input('Umidade relativa (%)......... '))
    q1=  float(input('Vazão de ar (m3/h)........... '))
    print('********* Ponto de Estado 2 **********')
    tbs2=float(input('Temperatura bulbo seco(ºC)... '))
    ur2= float(input('Umidade relativa (%)......... '))
    q2=  float(input('Vazão de ar (m3/h)........... '))
    ur1 = ur1/100.
    if ur1 == 1.:
        ur1 = 0.99999
    ur2 = ur2/100.
    if ur2 == 1.:
        ur2 = 0.99999
    # Ponto de Estado 1
    pvs1= pressao_vapor_saturado(tbs1)
    pv1=  ur1* pvs1
    rm1=  razao_mistura1(pv1)
    e1=   entalpia(tbs1,rm1)
    ve1=  volume_especifico(tbs1,rm1)
    tbm1 = temperatura_b_molhado(tbs1, e1)
    tpo1 = temperatura_ponto_orvalho(pv1)
    #   Output  Ponto de Estado
    urr1=ur1*100.
    rma1=rm1*1000.
    print('************** Ponto de Estado 1 *********************')
    print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs1)
    print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm1)
    print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo1)
    print(' Umidade relativa(%%)..................... %7.3f'% urr1)
    print(' Razão de mistura (g/kg) ................ %7.3f' % rma1)
    print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv1)
    print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs1)
    print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
    print(' Entalpia(kJ/kg) ........................ %7.3f' % e1)
    print(' Volume específico(m3/kg) ............... %7.3f' % ve1)
    print('  ***** Vazão (m3/h) .................... %7.3f' % q1)
    # Ponto de Estado 2
    pvs2=pressao_vapor_saturado(tbs2)
    pv2= ur2*pvs2
    rm2= razao_mistura1(pv2)
    e2=entalpia(tbs2,rm2)
    ve2=volume_especifico(tbs2,rm2)
    tbm2=temperatura_b_molhado(tbs2,e2)
    tpo2=temperatura_ponto_orvalho(pv2)
    # Output - Ponto de Estado 2
    urr2=ur2*100.
    rma2=rm2*1000.
    print('\n************** Ponto de Estado 2 *******************')
    print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs2)
    print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm2)
    print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo2)
    print(' Umidade relativa(%%)..................... %7.3f'% urr2)
    print(' Razão de mistura (g/kg) ................ %7.3f' % rma2)
    print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv2)
    print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs2)
    print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
    print(' Entalpia(kJ/kg) ........................ %7.3f' % e2)
    print(' Volume específico(m3/kg) ............... %7.3f' % ve2)
    print('  ***** Vazão (m3/h) .................... %7.3f' % q2)
    # Ponto de Estado 3 -  Mistura de dois fluxos de ar
    q1 = q1 * ve1
    q2 = q2 * ve2
    tbs3=(q1*tbs1+q2*tbs2)/(q1+q2)
    rm3= (q1*rm1+q2*rm2)/(q1+q2)
    e3=  (q1*e1+q2*e2)  /(q1+q2)
    eref=e3
    pvs3=pressao_vapor_saturado(tbs3)
    pv3=pressao_vapor(rm3)
    ur3=pv3/pvs3
    if ur3 >= 1.00:
        print('Atenção! Formação de Neblina')
        ur3=0.99999
        while True:
            tbs3=tbs3+0.01
            pvs3=pressao_vapor_saturado(tbs3)
            rm3=razao_mistura1(pvs3)
            e3=entalpia(tbs3,rm3)
            if abs(eref-e3) == 0.0001:
                break
    else:
        ve3=volume_especifico(tbs3,rm3)
        tpo3=temperatura_ponto_orvalho(pv3)
        tbm3=temperatura_b_molhado(tbs3,e3)
        q3=(q1+q2)/ve3
        # Output - Ponto de Estado da Nistura
        ur3=ur3*100.
        rm3=rm3*1000.
        print('\n************** Ponto de Estado Mistura ************')
        print(' Temperatura de bulbo seco(ºC) .......... %7.3f' % tbs3)
        print(' Temperatura de bulbo molhado(ºC) ....... %7.3f' % tbm3)
        print(' Temperatura do ponto de orvalho(ºC) .... %7.3f' % tpo3)
        print(' Umidade relativa(%%)..................... %7.3f' % ur3)
        print(' Razão de mistura (g/kg) ................ %7.3f' % rm3)
        print(' Pressão parcial de vapor(kPa)........... %7.3f' % pv3)
        print(' Pressão de vapor de saturação(kPa) ..... %7.3f' % pvs3)
        print(' Pressão barométrica(kPa) ............... %7.3f' % patm)
        print(' Entalpia(kJ/kg) ........................ %7.3f' % e3)
        print(' Volume específico(m3/kg) ............... %7.3f' % ve3)
        print('  ***** Vazão total(m3/h) ............... %7.3f' % q3)

def qual_ponto():
    # cálculo Ponto de Estado  - input: tbs com outra propriedade
    global tbs, tbm, patm, tbs, ur, pvs, pv, rm, tpo, ve
    global e
    print('\n Ponto de Estado - Dados conhecidos: ')
    print('                  1. tbs e ur    ')
    print('                  2. tbs e tbm   ')
    print('                  3. tbs e tpo   ')
    par=int(input('Qual sua opção? '))
    if par == 1:
        pe_tbs_ur()
    elif par == 2:
        pe_tbs_tbm()
    else:
        pe_tbs_tpo()
    # Output
    ur=ur*100.
    rm=rm*1000.
    print('')
    print('*************** Ponto de Estado ***************')
    print('Temperatura de bulbo seco (ºC)......... %7.3f' %tbs)
    print('Temperatura de bulbo molhado (ºC)...... %7.3f' %tbm)
    print('Temperatura de ponto de orvalho (ºC) .. %7.3f' %tpo)
    print('Umidade relativa (%%) ................. %8.3f' %ur)
    print('Razão de mistura (g/kg)................ %7.3f' %rm)
    print('Pressão barométrica (kPa).............. %7.3f' %patm)
    print('Pressão de vapor saturado (kPa)........ %7.3f' %pvs)
    print('Pressão parcial de vapor (kPa)......... %7.3f' %pv)
    print('Entalpia (kJ/kg)....................... %7.3f' %e)
    print('Volume específico (m3/kg).............. %7.3f' %ve)

def qual_processo():
    global proc
    print('\n************************* Processos   ******************************* ')
    print('               1. aquecimento/resfriamento  ')
    print('               2. umificação adiabática:dado temperatura final')
    print('               3. umificação adiabática:dado umidade relativa final')
    print('               4. umificação adiabática:dado razão de mistura final')
    print('               5. mistura de dois fluxos de ar')
    print('********************************************************************* ')
    proc=int(input('Qual processo? '))
    if proc == 1:
        aquece_resfria()
    elif proc == 2:
        u_adiabatica_tbs()
    elif proc == 3:
        u_adiabatica_ur()
    elif proc == 4:
        u_adiabatica_rm()
    else:
        mistura_fluxos()
#
# ############################################# fim definição de funçôes #######
#
# ###################### programa principal  ###################################
if __name__ == '__main__':
    import numpy as np
    # apresentação da versão
    print_hi('\n Programa GRAPSI - versão 1.0 python \n')
    # calcula a pressão atmosférica em função da altitude
    patm = p_atm()
    # tipo de cálculo
    calc=int(input('\n Calcular===>  1.Ponto de estado ou 2.Processos?   '))
    if calc == 1:
        qual_ponto()
    if calc == 2:
        qual_processo()
