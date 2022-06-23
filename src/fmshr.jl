@fenicsclass Geometry

#functions necessary for creating meshes from geometrical objects.
#2D objects below
Circle(centre, radius) = Geometry(mshr.Circle(centre, radius))
Rectangle(corner1, corner2) = Geometry(mshr.Rectangle(corner1, corner2))
function Ellipse(centre, horizontal_semi_axis, vertical_semi_axis, fragments)
    Geometry(mshr.Ellipse(centre, horizontal_semi_axis, vertical_semi_axis, fragments))
end

#3d objects below
Box(corner1, corner2) = Geometry(mshr.Box(corner1, corner2))
function Cone(top, bottom, bottom_radius, slices::Int)
    Geometry(mshr.Cone(top, bottom, bottom_radius, slices))
end
Sphere(centre, radius) = Geometry(mshr.Sphere(centre, radius))

#we generate the mesh of some geometrical objects
function generate_mesh(geom_object::Geometry, size::Int)
    Mesh(mshr.generate_mesh(geom_object.pyobject, size))
end
export Circle, Rectangle, Ellipse, Box, Cone, Sphere, generate_mesh

#operator overloading for Geometry types so we can create composite shapes
function +(geom_object1::Geometry, geom_object2::Geometry)
    Geometry(geom_object1.pyobject.__add__(geom_object2.pyobject))
end
function -(geom_object1::Geometry, geom_object2::Geometry)
    Geometry(geom_object1.pyobject.__sub__(geom_object2.pyobject))
end
function *(geom_object1::Geometry, geom_object2::Geometry)
    Geometry(geom_object1.pyobject.__mul__(geom_object2.pyobject))
end

function set_subdomain(object::Geometry, degree, domain2)
    fenicspycall(object, :set_subdomain, degree, domain2.pyobject)
end
"""
Extrude2D

Extrudes a 2D geometry to 3D
"""
function Extrude2D(object::Geometry, thickness)
    Geometry(mshr.Extrude2D(object.pyobject, thickness))
end
export set_subdomain, Extrude2D

function CSGUnion(geom_object1::Geometry, geom_object2::Geometry)
    Geometry(mshr.CSGUnion(geom_object1.pyobject, geom_object2.pyobject))
end
function CSGIntersection(geom_object1::Geometry, geom_object2::Geometry)
    Geometry(mshr.CSGIntersection(geom_object1.pyobject, geom_object2.pyobject))
end
function CSGDifference(geom_object1::Geometry, geom_object2::Geometry)
    Geometry(mshr.CSGDifference(geom_object1.pyobject, geom_object2.pyobject))
end

function CSGScaling(geom_object1::Geometry, scale)
    Geometry(mshr.CSGScaling(geom_object1.pyobject, scale))
end

function CSGRotation(geom_object1::Geometry, angle)
    Geometry(mshr.CSGRotation(geom_object1.pyobject, angle))
end
function CSGRotation(geom_object1::Geometry, point, angle)
    Geometry(mshr.CSGRotation(geom_object1.pyobject, point, angle))
end
function CSGRotation(geom_object1::Geometry, point_axis, point_center, angle)
    Geometry(mshr.CSGRotation(geom_object1.pyobject, point_axis, point_center, angle))
end

function CSGTranslation(geom_object1::Geometry, translation_point)
    Geometry(mshr.CSGTranslation(geom_object1.pyobject, translation_point))
end

export CSGUnion, CSGIntersection, CSGDifference, CSGScaling, CSGRotation, CSGTranslation
