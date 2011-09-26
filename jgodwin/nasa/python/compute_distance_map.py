import vtk,numpy,triangleio,time

args = triangleio.initArguments()

outputname = args.get('out',None,required=True)
dx = args.get('dx',0.05,False,float)
filename = args.get('in',None,required=True)

pointsource = vtk.vtkProgrammableSource()

def readPoints():
    nv,vertices = triangleio.read(filename)
        
    output = pointsource.GetPolyDataOutput()
    points = vtk.vtkPoints()
    output.SetPoints(points)

    for vertex in vertices:
        x,y,z = vertex[0],vertex[1],vertex[2]
        points.InsertNextPoint(x,y,z)

pointsource.SetExecuteMethod(readPoints)

surf = vtk.vtkSurfaceReconstructionFilter()
surf.SetInputConnection(pointsource.GetOutputPort())
surf.SetSampleSpacing(dx)

writer = vtk.vtkStructuredPointsWriter()
writer.SetInputConnection(surf.GetOutputPort())
writer.SetFileTypeToASCII()
writer.SetFileName(outputname)
writer.Write()
