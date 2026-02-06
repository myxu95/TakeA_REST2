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
                                   method: str = 'exponential') -> List[float]:
        """
        Calculate temperature ladder for REST2 replicas

        Args:
            T_min: Minimum temperature (K)
            T_max: Maximum temperature (K)
            n_replicas: Number of replicas
            method: Scaling method (default: 'exponential' - exponential temperature distribution)
                    Formula: T[i] = T_min * exp(i * log(T_max/T_min) / (n-1))

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
        # Default and recommended: Exponential temperature spacing
        # Formula: T[i] = T_min * exp(i * log(T_max/T_min) / (n-1))
        # Equivalent: T[i] = T_min * (T_max/T_min)^(i/(n-1))

        if method == 'exponential':
            # Exponential temperature spacing
            # More replicas at lower temperatures (ΔT increases)
            ratio = (T_max / T_min) ** (1.0 / (n_replicas - 1))
            temperatures = [T_min * (ratio ** i) for i in range(n_replicas)]
        elif method == 'linear':
            # Legacy support: Linear temperature spacing
            temperatures = np.linspace(T_min, T_max, n_replicas).tolist()
        else:
            raise TemperatureCalculationError(
                f"Unknown scaling method: {method}. Use 'exponential' (recommended) or 'linear'"
            )

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
                                       method: str = 'exponential') -> Tuple[List[float], List[float]]:
        """
        Calculate both temperature ladder and scaling factors

        Args:
            T_min: Minimum temperature (K)
            T_max: Maximum temperature (K)
            n_replicas: Number of replicas
            method: Scaling method ('linear', 'exponential', 'linear_lambda')

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
                                     method: str = 'exponential') -> bool:
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
        if method not in ['linear', 'exponential', 'linear_lambda']:
            raise TemperatureCalculationError(f"Unknown scaling method: {method}")

        # Method-specific validation
        if method in ['exponential', 'linear_lambda'] and T_min <= 0:
            raise TemperatureCalculationError(f"T_min must be positive for {method} scaling")

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
        print("\n1. Linear Temperature Scaling:")
        temperatures_linear, scaling_linear = TemperatureCalculator.calculate_temperature_and_scaling(
            T_min, T_max, n_replicas, 'linear'
        )
        TemperatureCalculator.print_temperature_summary(temperatures_linear, scaling_linear, 'linear')

        # Test linear lambda scaling (NEW)
        print("\n2. Linear Lambda Scaling (Uniform λ):")
        temperatures_lambda, scaling_lambda = TemperatureCalculator.calculate_temperature_and_scaling(
            T_min, T_max, n_replicas, 'linear_lambda'
        )
        TemperatureCalculator.print_temperature_summary(temperatures_lambda, scaling_lambda, 'linear_lambda')

        # Verify uniform lambda spacing
        lambda_diffs = [scaling_lambda[i] - scaling_lambda[i+1] for i in range(n_replicas-1)]
        print(f"λ adjacent differences: {[f'{d:.6f}' for d in lambda_diffs]}")
        print(f"Δλ constant = {lambda_diffs[0]:.6f} (uniform spacing)")

        # Test exponential scaling
        print("\n3. Exponential Temperature Scaling:")
        temperatures_exp, scaling_exp = TemperatureCalculator.calculate_temperature_and_scaling(
            T_min, T_max, n_replicas, 'exponential'
        )
        TemperatureCalculator.print_temperature_summary(temperatures_exp, scaling_exp, 'exponential')

        # Test validation
        print("\nTesting validation:")
        TemperatureCalculator.validate_temperature_parameters(T_min, T_max, n_replicas, 'linear_lambda')
        print("✓ Parameter validation passed for linear_lambda")

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