function makeVideoDemo()
    importer = LibraryImporter();
    importer.import('MeshProcessing')
    importer.import('StatisticalShapeModel')
    importer.import('Tools')
    importer.addFolder('..')

    template = shiftMeshToOrigin(loadSTL('template.stl'));
    template = transformMesh(template, [7 7 7], rotationMatrix('z', 10/180*pi));
    target = shiftMeshToOrigin(loadSTL('target.stl'));

    makeVideo('demo.avi', template, target, ...
              'View', [4 6], ...
              'alpha', [1e9 1e5 1e3 100 10 1 0.1 0.01 0.001]')
end
