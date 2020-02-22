from abc import ABC, abstractmethod

class system(ABC):

    @abstractmethod
    def getAcceleration(self, particle):
        pass

    @abstractmethod
    def update(self, deltaT):
        pass