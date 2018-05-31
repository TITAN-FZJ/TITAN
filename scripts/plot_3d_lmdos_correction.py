#!/usr/bin/env python
import   numpy               as  np
import   matplotlib.pyplot   as  plt
import   matplotlib          as  mpl
from     matplotlib.patches import Rectangle

#####################################################
########### This is for the fonts and text ##########
#####################################################

mpl.rcParams['text.usetex']        = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.size']          = 14
mpl.rcParams['font.family']        = 'serif'
mpl.rcParams['font.serif']         = 'Computer Modern Roman'

########################################################################
############################# 3d LMDOS #################################
########################################################################
###### Cr ########
##################

file = open('New_ldos_BS_Bi2Te3/ldos0001.dat_Cr_Bulk')
lines_Cr=file.readlines()
Energy1=[]
ldos_dn_Cr =[]
ldos_up_Cr =[]

###############################################
############# Select needed data ##############
###############################################

for line in lines_Cr[3:]:
    data_str=line.split()
    Energy1.append(float(data_str[0]))
    ldos_dn_Cr.append(float(data_str[3]))
    ldos_up_Cr.append(float(data_str[5]))

################################################
############## Adjusting the units #############
################################################
Ef       = 0.5402747990
Energy1   = [(Energy1[i]-Ef)*13.60 for i in range(len(Energy1))]
ldos_dn_Cr  = [-ldos_dn_Cr[i]/13.60  for i in range(len(Energy1))]
ldos_up_Cr  = [ ldos_up_Cr[i]/13.60  for i in range(len(Energy1))]

##################
###### Mn ########
##################

file     = open('New_ldos_BS_Bi2Te3/ldos0001.dat_Cr_Surface')
lines_Mn = file.readlines()
Energy2=[]
ldos_dn_Mn =[]
ldos_up_Mn =[]

###############################################
############# Select needed data ##############
###############################################

for line in lines_Mn[3:]:
    data_str=line.split()
    Energy2.append(float(data_str[0]))
    ldos_dn_Mn.append(float(data_str[3]))
    ldos_up_Mn.append(float(data_str[5]))

################################################
############## Adjusting the units #############
################################################

Ef          = 0.6628189680
Energy2      = [(Energy2[i]-Ef)*13.60 for i in range(len(Energy2))]
ldos_dn_Mn  = [ ldos_dn_Mn[i]/13.60  for i in range(len(Energy2))]
ldos_up_Mn  = [-ldos_up_Mn[i]/13.60  for i in range(len(Energy2))]

############################################################################

fig, ax1 = plt.subplots()
ax1.annotate('Band gap', xy=(0.1, 5.0), xytext=(0.45, 6.0), arrowprops=dict(facecolor='black', shrink=0.05))
ax1.annotate('Bound states', xy=(0.05, 2.0), xytext=(0.4, 3.0), arrowprops=dict(facecolor='black', shrink=0.05))
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((-0.15,-8.0), 0.25, 15, facecolor="paleturquoise"))

########################################################################
############################# 3d LMDOS #################################
########################################################################
###### Pd ########
##################

file = open('New_ldos_BS_Bi2Te3/ldos0001.dat_Pd_Bulk')
lines_Pd=file.readlines()
Energy1=[]
ldos_dn_Pd =[]
ldos_up_Pd =[]

###############################################
############# Select needed data ##############
###############################################

for line in lines_Pd[3:]:
    data_str=line.split()
    Energy1.append(float(data_str[0]))
    ldos_dn_Pd.append(float(data_str[3]))
    ldos_up_Pd.append(float(data_str[5]))

################################################
############## Adjusting the units #############
################################################

Ef          = 0.5402747990
Energy1     = [(Energy1[i]-Ef)*13.60 for i in range(len(Energy1))]
ldos_dn_Pd  = [-ldos_dn_Pd[i]/13.60  for i in range(len(Energy1))]
ldos_up_Pd  = [ ldos_up_Pd[i]/13.60  for i in range(len(Energy1))]

########################################################################
############################# 3d LMDOS #################################
########################################################################
###### Pd ########
##################

file = open('New_ldos_BS_Bi2Te3/ldos0001.dat_Pd_Surface')
lines_Pd_S=file.readlines()
Energy2=[]
ldos_dn_Pd_S =[]
ldos_up_Pd_S =[]
###############################################
############# Select needed data ##############
###############################################
for line in lines_Pd_S[3:]:
    data_str=line.split()
    Energy2.append(float(data_str[0]))
    ldos_dn_Pd_S.append(float(data_str[3]))
    ldos_up_Pd_S.append(float(data_str[5]))
################################################
############## Adjusting the units #############
################################################
Ef            = 0.6628189680
Energy2       = [(Energy2[i]-Ef)*13.60 for i in range(len(Energy2))]
ldos_dn_Pd_S  = [-ldos_dn_Pd_S[i]/13.60  for i in range(len(Energy2))]
ldos_up_Pd_S  = [ ldos_up_Pd_S[i]/13.60  for i in range(len(Energy2))]
############################################################################
############################ 3d ldos subplots ##############################
############################################################################
plt.plot(Energy1, ldos_up_Cr, '-',  color="black", label='Cr Bulk'   )
plt.plot(Energy1, ldos_dn_Cr, '--', color="black")
plt.plot(Energy1, ldos_up_Pd, '-',  color="blue", label='Pd Bulk'   )
plt.plot(Energy1, ldos_dn_Pd, '--', color="blue" )
plt.plot(Energy2, ldos_dn_Mn, '-', color="red"  , label='Cr Surface')
plt.plot(Energy2, ldos_up_Mn, '--', color="red"  )
plt.plot(Energy2, ldos_dn_Pd_S, '--',  color="Green"  )
plt.plot(Energy2, ldos_up_Pd_S, '-', color="Green" , label='Pd Bulk')
####################
#### Line at Ef ####
####################
plt.axvline(x=0, color='black')
plt.axhline(y=0, color='black')
plt.axis([-5.0, 2.0, -8.0, 7.0])
plt.xlabel('$\\varepsilon-\\varepsilon_\mathrm{F}$ (eV)')
plt.ylabel('LDOS (States/eV)')
plt.legend(loc=2, ncol=1, borderaxespad=0.0)
############################################################################
############################ Figures Format ################################
############################################################################
fig = mpl.pyplot.gcf()
fig.set_size_inches(7.0, 5.0)
plt.savefig('lmdos_Bulk_vs_surface.pdf')
plt.show()

############################################################################

