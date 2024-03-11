-- Lua (Keep this comment, this is an indication for editor's 'run' command)
S = scene_graph.find_or_create_object('OGF::MeshGrob', 'bunnin')
S.shader.painting  = 'ATTRIBUTE'
S.shader.attribute = 'facets.patch'
S.shader.autorange()