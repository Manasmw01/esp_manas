// Copyright (c) 2011-2024 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "kalman_filter_sysc_catapult.h"

#define DRV_NAME	"kalman_filter_sysc_catapult"

/* <<--regs-->> */
#define KALMAN_FILTER_MAC_N_REG 0x48
#define KALMAN_FILTER_MAC_VEC_REG 0x44
#define KALMAN_FILTER_MAC_LEN_REG 0x40

struct kalman_filter_sysc_catapult_device {
	struct esp_device esp;
};

static struct esp_driver kalman_filter_driver;

static struct of_device_id kalman_filter_device_ids[] = {
	{
		.name = "SLD_KALMAN_FILTER_SYSC_CATAPULT",
	},
	{
		.name = "eb_056",
	},
	{
		.compatible = "sld,kalman_filter_sysc_catapult",
	},
	{ },
};

static int kalman_filter_devs;

static inline struct kalman_filter_sysc_catapult_device *to_kalman_filter(struct esp_device *esp)
{
	return container_of(esp, struct kalman_filter_sysc_catapult_device, esp);
}

static void kalman_filter_prep_xfer(struct esp_device *esp, void *arg)
{
	struct kalman_filter_sysc_catapult_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->mac_n, esp->iomem + KALMAN_FILTER_MAC_N_REG);
	iowrite32be(a->mac_vec, esp->iomem + KALMAN_FILTER_MAC_VEC_REG);
	iowrite32be(a->mac_len, esp->iomem + KALMAN_FILTER_MAC_LEN_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool kalman_filter_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct kalman_filter_sysc_catapult_device *kalman_filter = to_kalman_filter(esp); */
	/* struct kalman_filter_sysc_catapult_access *a = arg; */

	return true;
}

static int kalman_filter_probe(struct platform_device *pdev)
{
	struct kalman_filter_sysc_catapult_device *kalman_filter;
	struct esp_device *esp;
	int rc;

	kalman_filter = kzalloc(sizeof(*kalman_filter), GFP_KERNEL);
	if (kalman_filter == NULL)
		return -ENOMEM;
	esp = &kalman_filter->esp;
	esp->module = THIS_MODULE;
	esp->number = kalman_filter_devs;
	esp->driver = &kalman_filter_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	kalman_filter_devs++;
	return 0;
 err:
	kfree(kalman_filter);
	return rc;
}

static int __exit kalman_filter_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct kalman_filter_sysc_catapult_device *kalman_filter = to_kalman_filter(esp);

	esp_device_unregister(esp);
	kfree(kalman_filter);
	return 0;
}

static struct esp_driver kalman_filter_driver = {
	.plat = {
		.probe		= kalman_filter_probe,
		.remove		= kalman_filter_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = kalman_filter_device_ids,
		},
	},
	.xfer_input_ok	= kalman_filter_xfer_input_ok,
	.prep_xfer	= kalman_filter_prep_xfer,
	.ioctl_cm	= KALMAN_FILTER_SYSC_CATAPULT_IOC_ACCESS,
	.arg_size	= sizeof(struct kalman_filter_sysc_catapult_access),
};

static int __init kalman_filter_init(void)
{
	return esp_driver_register(&kalman_filter_driver);
}

static void __exit kalman_filter_exit(void)
{
	esp_driver_unregister(&kalman_filter_driver);
}

module_init(kalman_filter_init)
module_exit(kalman_filter_exit)

MODULE_DEVICE_TABLE(of, kalman_filter_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("kalman_filter_sysc_catapult driver");
