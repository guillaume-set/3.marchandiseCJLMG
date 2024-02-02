import numpy as np
import matplotlib.pyplot
from itertools import combinations
from networkx import DiGraph
import requests
from vrpy import VehicleRoutingProblem


coordonnees =  {'dépot SORAGIS' : (2.36815, 48.74991), 'Mines Paris' : (2.33969, 48.84563),
'Observatoire de Paris' : (2.33650, 48.83730),
'Marie du 14e' : (2.32698, 48.83320),
'Gare Montparnasse TGV' : (2.32159, 48.84117),
'Mairie du 15e' : (2.29991, 48.84126)}

def matrice_adj (coords) :
    M = np.zeros((6,6))
    n = 0
    j = 0
    for key1 in coords :
        for key2 in coords :
            (lng_origine, lat_origine) = (coords[key1][0], coords[key1][1])
            (lng_destination, lat_destination) = (coords[key2][0], coords[key2][1])
            try:
                r = requests.get(f"https://wxs.ign.fr/essentiels/geoportail/itineraire/rest/1.0.0/route?resource=bdtopo-osrm&start={lng_origine},{lat_origine}&end={lng_destination},{lat_destination}").json()
            except Exception :
                print(f'erreur requete !')
            M[n][j] = (r['duration'])
            j = j+1
            if key2 == 'Mairie du 15e' :
                j = 0
        n = n+1
    return M

matrice_distances = matrice_adj(coordonnees)


X = np.array([coord[0] for coord in coordonnees.values()])
Y = np.array([coord[1] for coord in coordonnees.values()])
G3 = DiGraph()
for (i, (xval, yval)) in enumerate(zip(X, Y)):
  G3.add_edge("Source", i, cost = matrice_distances[0][i])
  G3.add_edge(i, "Sink", cost = matrice_distances[0][i])
  G3.nodes[i]["demand"] = 1
for ((i, (x1, y1)), (j, (x2, y2))) in combinations(list(enumerate(zip(X, Y))), 2):
  G3.add_edge(i, j, cost = matrice_distances[i][j])
  G3.add_edge(j, i, cost = matrice_distances[j][i])
probleme_optimisation = VehicleRoutingProblem(G3, load_capacity = 9999)
probleme_optimisation.solve()
matplotlib.pyplot.figure()
matplotlib.pyplot.scatter(X, Y, color = 'red', marker = 'x')
for route in probleme_optimisation.best_routes.values():
  itineraire_X =[X[p] for p in route[1:-1]] + [2.36815]
  itineraire_Y =[Y[p] for p in route[1:-1]] + [48.74991]
  for (a, b, c, d) in zip(itineraire_X[:-1], itineraire_Y[:-1], [x2-x1 for (x1, x2) in zip(itineraire_X[:-
1], itineraire_X[1:])], [y2-y1 for (y1, y2) in zip(itineraire_Y[:-1], itineraire_Y[1:])]):
    matplotlib.pyplot.arrow(a, b, c, d, width = 0.0003, color = 'black', length_includes_head =
True, head_width = 0.002)
matplotlib.pyplot.show()




#Calcul énergétique

#N1= 3.5 tonnes
#N2 = 5.2 tonnes
#N3= 18 tonnes
#N4 = 44 tonnes
#le dico caracteristiques a pour clé "modèle de voiture", et a comme valeur [masse, aire frontale, coeff Cx, coeff k1]
caracteristiques= {'N1':[3.5e3, 4.56, 0.46e-3,2.15e-3], 'N2':[5.2e3, 4.86, 0.53e-3,2.5e-3], 'N3':[18e3, 7.17, 0.69e-3,3.25e-3], 'N4':[44e3, 9.34, 0.79e-3,3.75e-3]}
rho=1.2 #kg/m^3
g= 9.81 #m/s^2
v= 9.7 #m/s = 35 km/h #on détermine une vitesse moyenne, tous les véhicules roulent à cette vitesse
M= matrice_adj(coordonnees)

destinations={'Mines Paris' : 1,
'Observatoire de Paris' : 2,
'Marie du 14e' : 3,
'Gare Montparnasse TGV' : 4,
'Mairie du 15e' : 5,'dépot SORAGIS': 0}

def pertes_nrj(depart,arrivee, vehicule):
    i= destinations[depart]
    j= destinations[arrivee]
    duree= M[i][j]
    m= caracteristiques[vehicule][0]
    A= caracteristiques[vehicule][1]
    Cx= caracteristiques[vehicule][2]
    k1= caracteristiques[vehicule][3]

    p=v*((rho*v**2*Cx*A)/2 + m*g*(k1 + k1*v**2)) #formule de la puissance
    E= p*duree #E est en Joules
    e=E/(3.6e+6)
    return e

import folium

m = folium.Map([48.866669,2.33333], zoom_start=12)

coordsls = {}
for key in coordonnees :
  coords_lieu = coordonnees[key]
  coordsls[key] = list(coords_lieu)

lieux = {}
cles = coordsls.keys()
ls_cles = list(cles)
for elem in ls_cles :
  lieux[f'{coordsls[elem]}'] = elem



itineraire=list(zip(itineraire_X,itineraire_Y))

E=0
for i in range(len(itineraire)-1):
  
    point=itineraire[i]
    point_2=itineraire[i+1]
    (lng_origine, lat_origine) = (point[0], point[1])
    (lng_destination, lat_destination) = (point_2[0], point_2[1])
    try:
            r = requests.get(f"https://wxs.ign.fr/essentiels/geoportail/itineraire/rest/1.0.0/route?resource=bdtopo-osrm&start={lng_origine},{lat_origine}&end={lng_destination},{lat_destination}").json()
    except Exception :
            print(f'erreur requete !')

    trail_coordinates = r['geometry']['coordinates']
    E=E + pertes_nrj(lieux[f'{list(point)}'],lieux[f'{list(point_2)}'],'N1' )
    for j in range(len(trail_coordinates)):
      k=trail_coordinates[j]
      x=(k[1],k[0])
      trail_coordinates[j]=x
    if i==0:
      folium.Marker(location=[point[1],point[0]],tooltip="Click me!",popup=lieux[f'{list(point)}'],icon=folium.Icon(color="red",icon="flag")).add_to(m)
    elif i!=0:
     folium.Marker(location=[point[1],point[0]],tooltip="Click me!",popup=lieux[f'{list(point)}'],icon=folium.Icon(color="green",icon="cloud")).add_to(m)

    folium.PolyLine(trail_coordinates, tooltip="Coast").add_to(m)

print(f'cela fera {E} kwh')
m.save('ma_carte.html')