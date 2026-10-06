"""Modal invariant-region bounds using exact relative-angle electrical ports."""
from modal_comparison_contract import Comparison
from graph_relative_model import GraphRelativeModel


class GraphComparison(Comparison):
    def __init__(self):
        super().__init__("dq")
        self.model=GraphRelativeModel(center=self.center["offset"])
    def evaluate(self,r,epsilon):
        self.model.set_support(self.Ti,self.groups,r)
        return super().evaluate(r,epsilon)
