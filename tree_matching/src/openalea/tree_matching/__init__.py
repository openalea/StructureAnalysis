from .__tree_matching__ import *

def build_treenode(id,father):
    return TreeNode(id,father)

TreeNode.factory().setBuilder(build_treenode)
del build_treenode