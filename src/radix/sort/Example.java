/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package radix.sort;

import java.util.Random;

public class Example {

	public static void main(String[] args) {

		Random rand = new Random(42L);
		for (int n = 0; n < 10; ++n) {
			int[] data = new int[32 * 1024 * 1024];
			for (int i = 0; i < data.length; ++i) {
				data[i] = rand.nextInt();
			}

			long t = System.currentTimeMillis();
			RadixSort.radixSort(data);
			//Arrays.sort(data);
			t = System.currentTimeMillis() - t;
			System.out.println("time: " + t + " ms");
			checkSorted(data);
		}
	}

	private static final void checkSorted(int[] data) {
		for (int i = 1; i < data.length; ++i) {
			if (data[i - 1] > data[i]) {
				throw new IllegalStateException("Not sorted!");
			}
		}
	}

}