#!/usr/bin/env python3
"""
REST2 Temperature Calculator Module
Unified temperature ladder and scaling factor calculations for REST2 simulations
"""

import numpy as np
from typing import List, Tuple, Union
import math


class TemperatureCalculationError(Exception):
    """Temperature calculation error"""
    pass


class TemperatureCalculator:
    """
    Unified temperature calculator for REST2 simulations
    Handles temperature ladder generation and scaling factor calculations
    """
    
    @staticmethod
    def calculate_temperature_ladder(T_min: float, T_max: float, n_replicas: int, 
                                   method: str = 'linear') -> List[float]:
        """
        Calculate temperature ladder for REST2 replicas
        
        Args:
            T_min: Minimum temperature (K)
            T_max: Maximum temperature (K)
            n_replicas: Number of replicas
            method: Scaling method ('linear' or 'exponential')
            
        Returns:
            List of temperatures for each replica
            
        Raises:
            TemperatureCalculationError: If parameters are invalid
        """
        # Validate input parameters
        if T_min <= 0:
            raise TemperatureCalculationError("T_min must be positive")
        if T_max <= T_min:
            raise TemperatureCalculationError("T_max must be greater than T_min")
        if n_replicas < 1:
            raise TemperatureCalculationError("n_replicas must be at least 1")
        
        # Single replica case
        if n_replicas == 1:
            return [T_min]
        
        # Calculate temperatures based on method
        if method == 'linear':
            temperatures = np.linspace(T_min, T_max, n_replicas).tolist()
        elif method == 'exponential':
            # Exponential temperature spacing
            ratio = (T_max / T_min) ** (1.0 / (n_replicas - 1))
            temperatures = [T_min * (ratio ** i) for i in range(n_replicas)]
        else:
            raise TemperatureCalculationError(f"Unknown scaling method: {method}")
        
        return temperatures
    
    @staticmethod
    def calculate_scaling_factors(temperatures: List[float]) -> List[float]:
        """
        Calculate REST2 scaling factors (λ = T_ref/T)
        
        Args:
            temperatures: List of temperatures for each replica
            
        Returns:
            List of scaling factors for each replica
            
        Raises:
            TemperatureCalculationError: If temperatures are invalid
        """
        if not temperatures:
            raise TemperatureCalculationError("Temperature list cannot be empty")
        
        if any(T <= 0 for T in temperatures):
            raise TemperatureCalculationError("All temperatures must be positive")
        
        # Reference temperature is the lowest temperature
        T_ref = min(temperatures)
        
        # Calculate scaling factors
        scaling_factors = [T_ref / T for T in temperatures]
        
        return scaling_factors
    
    @staticmethod
    def calculate_temperature_and_scaling(T_min: float, T_max: float, n_replicas: int,
                                       method: str = 'linear') -> Tuple[List[float], List[float]]:
        """
        Calculate both temperature ladder and scaling factors
        
        Args:
            T_min: Minimum temperature (K)
            T_max: Maximum temperature (K)
            n_replicas: Number of replicas
            method: Scaling method ('linear' or 'exponential')
            
        Returns:
            Tuple of (temperatures, scaling_factors)
        """
        temperatures = TemperatureCalculator.calculate_temperature_ladder(
            T_min, T_max, n_replicas, method
        )
        scaling_factors = TemperatureCalculator.calculate_scaling_factors(temperatures)
        
        return temperatures, scaling_factors
    
    @staticmethod
    def validate_temperature_parameters(T_min: float, T_max: float, n_replicas: int,
                                     method: str = 'linear') -> bool:
        """
        Validate temperature calculation parameters
        
        Args:
            T_min: Minimum temperature (K)
            T_max: Maximum temperature (K)
            n_replicas: Number of replicas
            method: Scaling method
            
        Returns:
            True if parameters are valid
            
        Raises:
            TemperatureCalculationError: If parameters are invalid
        """
        # Basic parameter validation
        if T_min <= 0:
            raise TemperatureCalculationError("T_min must be positive")
        if T_max <= T_min:
            raise TemperatureCalculationError("T_max must be greater than T_min")
        if n_replicas < 1:
            raise TemperatureCalculationError("n_replicas must be at least 1")
        if method not in ['linear', 'exponential']:
            raise TemperatureCalculationError(f"Unknown scaling method: {method}")
        
        # Method-specific validation
        if method == 'exponential' and T_min <= 0:
            raise TemperatureCalculationError("T_min must be positive for exponential scaling")
        
        return True
    
    @staticmethod
    def print_temperature_summary(temperatures: List[float], scaling_factors: List[float],
                                method: str = 'linear') -> None:
        """
        Print formatted temperature ladder summary
        
        Args:
            temperatures: List of temperatures
            scaling_factors: List of scaling factors
            method: Scaling method used
        """
        print(f"\nTemperature Ladder ({method} scaling):")
        print("-" * 50)
        print(f"{'Replica':<8} {'Temperature (K)':<15} {'Scaling Factor (λ)':<18}")
        print("-" * 50)
        
        for i, (T, lambda_val) in enumerate(zip(temperatures, scaling_factors)):
            print(f"{i:<8} {T:<15.1f} {lambda_val:<18.6f}")
        
        print("-" * 50)
        print(f"Temperature range: {min(temperatures):.1f} - {max(temperatures):.1f} K")
        print(f"Reference temperature: {min(temperatures):.1f} K")
        print()


def main():
    """Test temperature calculator"""
    try:
        # Test parameters
        T_min = 300.0
        T_max = 340.0
        n_replicas = 8
        
        print("Testing Temperature Calculator")
        print("=" * 40)
        
        # Test linear scaling
        print("Linear scaling:")
        temperatures_linear, scaling_linear = TemperatureCalculator.calculate_temperature_and_scaling(
            T_min, T_max, n_replicas, 'linear'
        )
        TemperatureCalculator.print_temperature_summary(temperatures_linear, scaling_linear, 'linear')
        
        # Test exponential scaling
        print("Exponential scaling:")
        temperatures_exp, scaling_exp = TemperatureCalculator.calculate_temperature_and_scaling(
            T_min, T_max, n_replicas, 'exponential'
        )
        TemperatureCalculator.print_temperature_summary(temperatures_exp, scaling_exp, 'exponential')
        
        # Test validation
        print("Testing validation:")
        TemperatureCalculator.validate_temperature_parameters(T_min, T_max, n_replicas, 'linear')
        print("✓ Parameter validation passed")
        
        # Test error handling
        print("\nTesting error handling:")
        try:
            TemperatureCalculator.calculate_temperature_ladder(0, T_max, n_replicas)
        except TemperatureCalculationError as e:
            print(f"✓ Caught expected error: {e}")
        
        print("\n✓ All tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main() 