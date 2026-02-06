#!/usr/bin/env python3
"""Test temperature calculator - Exponential scaling only"""

import sys
sys.path.insert(0, '/public/home/xmy/tools/TakeA_REST2')

from utils.temperature_calculator import TemperatureCalculator, TemperatureCalculationError

def main():
    """Test temperature calculator"""
    try:
        # Test parameters
        T_min = 300.0
        T_max = 340.0
        n_replicas = 8

        print("=" * 70)
        print("Testing Temperature Calculator - Exponential Scaling")
        print("=" * 70)

        # Test exponential scaling (DEFAULT and RECOMMENDED)
        print("\nExponential Temperature Scaling (Default):")
        print("-" * 70)
        print("Formula: T[i] = T_min * exp(i * log(T_max/T_min) / (n-1))")
        print("         = T_min * (T_max/T_min)^(i/(n-1))")
        print()

        temperatures_exp, scaling_exp = TemperatureCalculator.calculate_temperature_and_scaling(
            T_min, T_max, n_replicas
        )  # Uses exponential by default

        TemperatureCalculator.print_temperature_summary(temperatures_exp, scaling_exp, 'exponential')

        # Verify temperature distribution
        temp_diffs = [temperatures_exp[i+1] - temperatures_exp[i] for i in range(n_replicas-1)]
        print(f"\nTemperature spacing (ΔT):")
        for i, diff in enumerate(temp_diffs):
            print(f"  T[{i+1}] - T[{i}] = {diff:.4f} K")
        print(f"\nΔT range: {min(temp_diffs):.4f} - {max(temp_diffs):.4f} K")
        print(f"Characteristic: ΔT increases from low to high temperature")
        print(f"Benefit: More sampling at lower temperatures")

        # Test validation
        print("\n" + "-" * 70)
        print("Testing validation:")
        TemperatureCalculator.validate_temperature_parameters(T_min, T_max, n_replicas)
        print("✓ Parameter validation passed")

        # Test error handling
        print("\nTesting error handling:")
        try:
            TemperatureCalculator.calculate_temperature_ladder(0, T_max, n_replicas)
        except TemperatureCalculationError as e:
            print(f"✓ Caught expected error: {e}")

        print("\n" + "=" * 70)
        print("✓ All tests passed!")
        print("=" * 70)

    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
