-- Lua (Keep this comment, this is an indication for editor's 'run' command)
S = scene_graph.find_or_create_object('OGF::MeshGrob', 'output')
S.shader.painting  = 'ATTRIBUTE'
-- S.shader.attribute = 'facets.patch'
S.shader.mesh_style='true; 0 0 0 1; 1'
S.shader.vertices_style='true; 0 1 1 1; 2'
S.query_interface('OGF::MeshGrobAttributesCommands').compute_vertices_id()
S.shader.attribute = 'vertices.id'
S.shader.autorange()
