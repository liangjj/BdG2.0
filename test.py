from parameters import Material, Geometry, Parameters

sn = Material(148.0, 0.000618, 0.66666, 0.1)
geometry = Geometry(200, 200, 30)
params = Parameters(sn, geometry)
