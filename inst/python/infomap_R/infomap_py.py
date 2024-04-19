from infomap import Infomap
import numpy

def infomap_python(matrix, n_clusters):
  im = Infomap("--flow-model undirected --preferred-number-of-modules "+str(n_clusters),silent=True)
  im.add_links(matrix)
  im.run()
  res = []
  for node in im.tree:
    if node.is_leaf:
      res.append(node.module_id)
  return res


