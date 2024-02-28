import heapq


class PrioritySet(object):
    def __init__(self):
        """
        Initialize the object with an empty heap and an empty set.
        """
        self.heap = []
        self.set = set()

    def add(self, d, pri):
        """
        Add an element to the heap if it is not already in the set.

        Args:
            d: The element to be added to the heap.
            pri: The priority of the element.
        """
        if d not in self.set:
            heapq.heappush(self.heap, (pri, d))
            self.set.add(d)

    def get(self):
        """
        Get the top element from the heap, remove it from the set, and return it.
        """
        pri, d = heapq.heappop(self.heap)
        self.set.remove(d)
        return d

    def __len__(self):
        return len(self.heap)
