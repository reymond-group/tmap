from typing import List, Dict, Any, Optional, Union
import uuid

class SynthesisTreeNode:
    """
    A node in a synthesis route tree.
    Represents a chemical compound or reaction in a synthesis route.
    """
    def __init__(
        self, 
        id: Optional[str] = None,
        label: str = "",
        data: Dict[str, Any] = None,
        children: List["SynthesisTreeNode"] = None
    ):
        """
        Initialize a synthesis tree node.
        
        Args:
            id: Unique identifier for the node, generated if not provided
            label: Label for the node (e.g., compound name)
            data: Additional data associated with the node
            children: List of child nodes
        """
        self.id = id if id is not None else str(uuid.uuid4())
        self.label = label
        self.data = data if data is not None else {}
        self.children = children if children is not None else []
        
    def add_child(self, child: "SynthesisTreeNode") -> None:
        """
        Add a child node to this node.
        
        Args:
            child: The child node to add
        """
        self.children.append(child)
        
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the node and its subtree to a dictionary representation.
        
        Returns:
            Dict: Dictionary representation of the node and its subtree
        """
        return {
            "id": self.id,
            "label": self.label,
            "data": self.data,
            "children": [child.to_dict() for child in self.children]
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SynthesisTreeNode":
        """
        Create a node and its subtree from a dictionary representation.
        
        Args:
            data: Dictionary representation of the node and its subtree
            
        Returns:
            SynthesisTreeNode: The constructed node and its subtree
        """
        node = cls(
            id=data.get("id"),
            label=data.get("label", ""),
            data=data.get("data", {})
        )
        
        for child_data in data.get("children", []):
            node.add_child(cls.from_dict(child_data))
            
        return node
    
    def __repr__(self) -> str:
        return f"SynthesisTreeNode(id={self.id}, label={self.label}, children={len(self.children)})"


class SynthesisTree:
    """
    A tree representing a chemical synthesis route.
    The root node typically represents the target compound,
    and children represent precursors or reactions.
    """
    def __init__(self, root: Optional[SynthesisTreeNode] = None):
        """
        Initialize a synthesis route tree.
        
        Args:
            root: Root node of the tree
        """
        self.root = root if root is not None else SynthesisTreeNode()
        
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the tree to a dictionary representation.
        
        Returns:
            Dict: Dictionary representation of the tree
        """
        return self.root.to_dict() if self.root else {}
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SynthesisTree":
        """
        Create a tree from a dictionary representation.
        
        Args:
            data: Dictionary representation of the tree
            
        Returns:
            SynthesisTree: The constructed tree
        """
        return cls(SynthesisTreeNode.from_dict(data))
    
    def __repr__(self) -> str:
        return f"SynthesisTree(root={self.root})" 